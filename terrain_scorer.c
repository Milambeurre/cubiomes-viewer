
/*
 * terrain_scorer.c
 *
 * Scoring topographique pour Cubiomes 26.2-snapshot-8.
 *
 * Objectif :
 *   - centre large et relativement plat
 *   - couronne de relief plus élevée / accidentée
 *   - relief présent dans plusieurs directions mais pas totalement fermé
 *   - grande zone extérieure relativement plate
 *   - petit bonus pour un sommet remarquable
 *
 * Le programme lit un fichier contenant une seed par ligne et produit un CSV
 * trié du meilleur au moins bon score.
 *
 * Compilation : voir README_FR.txt
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "generator.h"
#include "util.h"

#define STEP_BLOCKS       32
#define MAX_RADIUS_BLOCKS 1536
#define CENTER_RADIUS     320
#define RING_INNER        320
#define RING_OUTER        800
#define OUTER_INNER       800
#define OUTER_OUTER       1536

#define RING_SECTORS 8
#define OUTER_SECTORS 16

typedef struct {
    double sum;
    double sumsq;
    double minv;
    double maxv;
    long n;
} Stats;

typedef struct {
    uint64_t seed;
    double score;
    double center_flat;
    double ring_relief;
    double enclosure;
    double outer_flat;
    double contrast;
    double peak_bonus;

    double center_mean;
    double center_sd;
    double ring_mean;
    double ring_sd;
    double outer_mean;
    double outer_sd;

    int ring_good_sectors;
    int outer_best_run;
    double max_height;
} Result;

static void stats_init(Stats *s)
{
    s->sum = 0.0;
    s->sumsq = 0.0;
    s->minv = 1e30;
    s->maxv = -1e30;
    s->n = 0;
}

static void stats_add(Stats *s, double v)
{
    s->sum += v;
    s->sumsq += v * v;
    if (v < s->minv) s->minv = v;
    if (v > s->maxv) s->maxv = v;
    s->n++;
}

static double stats_mean(const Stats *s)
{
    return s->n ? s->sum / (double)s->n : 0.0;
}

static double stats_sd(const Stats *s)
{
    if (s->n < 2) return 0.0;
    double m = stats_mean(s);
    double v = s->sumsq / (double)s->n - m * m;
    return v > 0.0 ? sqrt(v) : 0.0;
}

static double clamp01(double x)
{
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

/* Score 0..100 : 100 <= good, 0 >= bad */
static double score_low(double x, double good, double bad)
{
    if (x <= good) return 100.0;
    if (x >= bad) return 0.0;
    return 100.0 * (1.0 - (x - good) / (bad - good));
}

/* Score 0..100 : 0 <= bad, 100 >= good */
static double score_high(double x, double bad, double good)
{
    if (x <= bad) return 0.0;
    if (x >= good) return 100.0;
    return 100.0 * (x - bad) / (good - bad);
}

/* Score en cloche : maximum dans [good_lo, good_hi]. */
static double score_band(double x, double bad_lo, double good_lo,
                         double good_hi, double bad_hi)
{
    if (x <= bad_lo || x >= bad_hi) return 0.0;
    if (x >= good_lo && x <= good_hi) return 100.0;
    if (x < good_lo)
        return 100.0 * (x - bad_lo) / (good_lo - bad_lo);
    return 100.0 * (bad_hi - x) / (bad_hi - good_hi);
}

static int sector_index(double angle, int n)
{
    /* atan2 => [-pi,pi], converti en [0,2pi). */
    double a = angle;
    if (a < 0.0) a += 2.0 * M_PI;
    int idx = (int)floor(a / (2.0 * M_PI) * n);
    if (idx < 0) idx = 0;
    if (idx >= n) idx = n - 1;
    return idx;
}

static int best_circular_run(const int *good, int n)
{
    int best = 0, cur = 0;
    for (int i = 0; i < 2 * n; i++) {
        if (good[i % n]) {
            cur++;
            if (cur > best) best = cur;
            if (best >= n) return n;
        } else {
            cur = 0;
        }
    }
    if (best > n) best = n;
    return best;
}

static double sample_height(Generator *g, SurfaceNoise *sn, int bx, int bz)
{
    /* mapApproxHeight travaille à l'échelle horizontale 1:4. */
    int qx = bx / 4;
    int qz = bz / 4;
    float y = 0.0f;
    if (mapApproxHeight(&y, NULL, g, sn, qx, qz, 1, 1) != 0)
        return NAN;
    return (double)y;
}

static Result score_seed(uint64_t seed)
{
    Result r;
    memset(&r, 0, sizeof(r));
    r.seed = seed;

    Generator g;
    setupGenerator(&g, MC_26_2_S8, 0);
    applySeed(&g, DIM_OVERWORLD, seed);

    SurfaceNoise sn;
    initSurfaceNoise(&sn, DIM_OVERWORLD, seed);

    Stats center, ring, outer;
    Stats ring_sec[RING_SECTORS];
    Stats outer_sec[OUTER_SECTORS];

    stats_init(&center);
    stats_init(&ring);
    stats_init(&outer);
    for (int i = 0; i < RING_SECTORS; i++) stats_init(&ring_sec[i]);
    for (int i = 0; i < OUTER_SECTORS; i++) stats_init(&outer_sec[i]);

    double max_h = -1e30;

    for (int z = -MAX_RADIUS_BLOCKS; z <= MAX_RADIUS_BLOCKS; z += STEP_BLOCKS) {
        for (int x = -MAX_RADIUS_BLOCKS; x <= MAX_RADIUS_BLOCKS; x += STEP_BLOCKS) {
            double d = sqrt((double)x * x + (double)z * z);
            if (d > MAX_RADIUS_BLOCKS) continue;

            double h = sample_height(&g, &sn, x, z);
            if (!isfinite(h)) continue;
            if (h > max_h) max_h = h;

            if (d <= CENTER_RADIUS) {
                stats_add(&center, h);
            } else if (d > RING_INNER && d <= RING_OUTER) {
                stats_add(&ring, h);
                int s = sector_index(atan2((double)z, (double)x), RING_SECTORS);
                stats_add(&ring_sec[s], h);
            } else if (d > OUTER_INNER && d <= OUTER_OUTER) {
                stats_add(&outer, h);
                int s = sector_index(atan2((double)z, (double)x), OUTER_SECTORS);
                stats_add(&outer_sec[s], h);
            }
        }
    }

    double cm = stats_mean(&center);
    double cs = stats_sd(&center);
    double rm = stats_mean(&ring);
    double rs = stats_sd(&ring);
    double om = stats_mean(&outer);
    double os = stats_sd(&outer);

    r.center_mean = cm;
    r.center_sd = cs;
    r.ring_mean = rm;
    r.ring_sd = rs;
    r.outer_mean = om;
    r.outer_sd = os;
    r.max_height = max_h;

    /*
     * 1) Centre plat.
     * Très bon <= 7 blocs d'écart-type ; très mauvais >= 24.
     */
    r.center_flat = score_low(cs, 7.0, 24.0);

    /*
     * 2) Contraste vertical couronne -> centre.
     * 0 si <= 6 blocs, 100 si >= 30.
     */
    double delta = rm - cm;
    r.contrast = score_high(delta, 6.0, 30.0);

    /*
     * 3) Relief de la couronne.
     * Combine dispersion et amplitude globale.
     */
    double ring_amp = ring.maxv - ring.minv;
    double relief_sd = score_high(rs, 8.0, 24.0);
    double relief_amp = score_high(ring_amp, 35.0, 95.0);
    r.ring_relief = 0.55 * relief_sd + 0.45 * relief_amp;

    /*
     * 4) "Amphithéâtre" :
     * un secteur est marqué si sa moyenne est >= centre + 12
     * OU si son max est >= centre + 35.
     *
     * On préfère 5 à 7 secteurs sur 8. 8/8 reste acceptable mais moins idéal
     * car cela ressemble davantage à une cuvette totalement fermée.
     */
    int good_ring = 0;
    for (int i = 0; i < RING_SECTORS; i++) {
        double sm = stats_mean(&ring_sec[i]);
        double sx = ring_sec[i].maxv;
        if (ring_sec[i].n > 0 &&
            ((sm - cm) >= 12.0 || (sx - cm) >= 35.0))
            good_ring++;
    }
    r.ring_good_sectors = good_ring;

    switch (good_ring) {
        case 5: case 6: case 7: r.enclosure = 100.0; break;
        case 4: r.enclosure = 78.0; break;
        case 8: r.enclosure = 72.0; break;
        case 3: r.enclosure = 55.0; break;
        case 2: r.enclosure = 30.0; break;
        case 1: r.enclosure = 12.0; break;
        default: r.enclosure = 0.0; break;
    }

    /*
     * 5) Extérieur constructible.
     * On ne demande PAS que toute la couronne 800-1536 soit plate.
     * On cherche un large arc continu de secteurs relativement plats.
     */
    int flat_outer[OUTER_SECTORS];
    for (int i = 0; i < OUTER_SECTORS; i++) {
        double sd = stats_sd(&outer_sec[i]);
        double amp = outer_sec[i].maxv - outer_sec[i].minv;

        /* Seuils volontairement souples au premier passage. */
        flat_outer[i] =
            (outer_sec[i].n > 0 && sd <= 14.0 && amp <= 55.0) ? 1 : 0;
    }
    int run = best_circular_run(flat_outer, OUTER_SECTORS);
    r.outer_best_run = run;

    /*
     * 16 secteurs => 22.5° chacun.
     * 5 secteurs contigus = 112.5° : excellent potentiel pour une grande zone.
     */
    if (run >= 5) r.outer_flat = 100.0;
    else if (run == 4) r.outer_flat = 82.0;
    else if (run == 3) r.outer_flat = 62.0;
    else if (run == 2) r.outer_flat = 38.0;
    else if (run == 1) r.outer_flat = 16.0;
    else r.outer_flat = 0.0;

    /*
     * 6) Sommet remarquable : bonus faible uniquement.
     */
    double peak_delta = max_h - cm;
    r.peak_bonus = score_high(peak_delta, 45.0, 95.0);

    /*
     * Pondération finale.
     * Les quatre critères structurants dominent ; le sommet ne peut jamais
     * sauver une mauvaise géographie.
     */
    r.score =
          0.24 * r.center_flat
        + 0.20 * r.contrast
        + 0.19 * r.ring_relief
        + 0.18 * r.enclosure
        + 0.16 * r.outer_flat
        + 0.03 * r.peak_bonus;

    return r;
}

static int cmp_result_desc(const void *a, const void *b)
{
    const Result *ra = (const Result *)a;
    const Result *rb = (const Result *)b;
    if (ra->score < rb->score) return 1;
    if (ra->score > rb->score) return -1;
    if (ra->seed < rb->seed) return -1;
    if (ra->seed > rb->seed) return 1;
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr,
            "Usage: %s seeds.txt classement.csv\n", argv[0]);
        return 1;
    }

    FILE *f = fopen(argv[1], "r");
    if (!f) {
        perror("Impossible d'ouvrir le fichier de seeds");
        return 1;
    }

    size_t cap = 2048, n = 0;
    uint64_t *seeds = (uint64_t *)malloc(cap * sizeof(uint64_t));
    if (!seeds) {
        fclose(f);
        fprintf(stderr, "Memoire insuffisante.\n");
        return 1;
    }

    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *p = line;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (*p == '\0' || *p == '#') continue;

        char *end = NULL;
        unsigned long long v = strtoull(p, &end, 10);
        if (end == p) continue;

        if (n == cap) {
            cap *= 2;
            uint64_t *tmp = (uint64_t *)realloc(seeds, cap * sizeof(uint64_t));
            if (!tmp) {
                free(seeds);
                fclose(f);
                fprintf(stderr, "Memoire insuffisante.\n");
                return 1;
            }
            seeds = tmp;
        }
        seeds[n++] = (uint64_t)v;
    }
    fclose(f);

    if (n == 0) {
        free(seeds);
        fprintf(stderr, "Aucune seed valide trouvee.\n");
        return 1;
    }

    Result *res = (Result *)malloc(n * sizeof(Result));
    if (!res) {
        free(seeds);
        fprintf(stderr, "Memoire insuffisante.\n");
        return 1;
    }

    fprintf(stderr, "Scoring de %zu seeds...\n", n);

    for (size_t i = 0; i < n; i++) {
        res[i] = score_seed(seeds[i]);

        if ((i + 1) % 25 == 0 || i + 1 == n) {
            fprintf(stderr, "\r%zu / %zu", i + 1, n);
            fflush(stderr);
        }
    }
    fprintf(stderr, "\nTri...\n");

    qsort(res, n, sizeof(Result), cmp_result_desc);

    FILE *o = fopen(argv[2], "w");
    if (!o) {
        perror("Impossible de creer le CSV");
        free(res);
        free(seeds);
        return 1;
    }

    fprintf(o,
        "rank,seed,score,center_flat,contrast,ring_relief,enclosure,"
        "outer_flat,peak_bonus,center_mean,center_sd,ring_mean,ring_sd,"
        "outer_mean,outer_sd,ring_good_sectors,outer_best_run,max_height\n");

    for (size_t i = 0; i < n; i++) {
        const Result *q = &res[i];
        fprintf(o,
            "%zu,%llu,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,"
            "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d,%.3f\n",
            i + 1,
            (unsigned long long)q->seed,
            q->score,
            q->center_flat,
            q->contrast,
            q->ring_relief,
            q->enclosure,
            q->outer_flat,
            q->peak_bonus,
            q->center_mean,
            q->center_sd,
            q->ring_mean,
            q->ring_sd,
            q->outer_mean,
            q->outer_sd,
            q->ring_good_sectors,
            q->outer_best_run,
            q->max_height);
    }

    fclose(o);

    fprintf(stderr, "Termine : %s\n", argv[2]);
    fprintf(stderr, "Top 10 :\n");
    for (size_t i = 0; i < n && i < 10; i++) {
        fprintf(stderr, "%2zu. seed=%llu score=%.2f\n",
            i + 1,
            (unsigned long long)res[i].seed,
            res[i].score);
    }

    free(res);
    free(seeds);
    return 0;
}
