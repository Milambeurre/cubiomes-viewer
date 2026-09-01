/*
 * terrain_scorer_v2.c
 *
 * Scorer topographique V2 pour Cubiomes 26.2-snapshot-8.
 *
 * Cible :
 *   - bassin central compact, ~250-350 blocs de diametre (r ~= 125-175)
 *   - centre regulier dans toutes les directions
 *   - debut du relief proche du centre
 *   - couronne principale de relief ~200-450 blocs
 *   - enceinte forte mais imparfaite (idealement 6-7 secteurs sur 8)
 *   - ouvertures peu nombreuses et contigues
 *   - espace terrestre substantiel au-dela du relief, ~450-1200 blocs
 *   - sommet remarquable = bonus uniquement
 *
 * Important : mapApproxHeight() est une approximation. Les meilleurs candidats
 * doivent etre controles visuellement puis dans Minecraft 26.2.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

#include "generator.h"
#include "util.h"
#include "biomes.h"

#define STEP_BLOCKS 24
#define MAX_RADIUS 1200

#define CORE_RADIUS 150
#define CORE_EDGE_INNER 120
#define CORE_EDGE_OUTER 175

#define TRANS_INNER 150
#define TRANS_OUTER 250

#define RING_INNER 200
#define RING_OUTER 450

#define OUTER_INNER 450
#define OUTER_OUTER 1200

#define SECTORS 8
#define OUTER_SECTORS 16

typedef struct {
    double sum, sumsq, minv, maxv;
    long n;
} Stats;

typedef struct {
    uint64_t seed;
    double score;

    double core_flat;
    double core_circular;
    double close_slope;
    double ring_strength;
    double enclosure;
    double openings;
    double outer_land;
    double peak_bonus;

    double core_mean, core_sd, core_range;
    double transition_mean, transition_rise;
    double ring_mean, ring_sd;
    int ring_good_sectors;
    int ring_best_good_run;
    int opening_groups;
    int largest_opening_run;

    double outer_land_fraction;
    int outer_land_sectors;
    int outer_best_land_run;
    double max_ring_height;
    double max_ring_rise;
} Result;

static void stats_init(Stats *s)
{
    s->sum = s->sumsq = 0.0;
    s->minv = 1e30;
    s->maxv = -1e30;
    s->n = 0;
}
static void stats_add(Stats *s, double v)
{
    s->sum += v; s->sumsq += v*v;
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
    double v = s->sumsq / (double)s->n - m*m;
    return v > 0.0 ? sqrt(v) : 0.0;
}
static double score_low(double x, double good, double bad)
{
    if (x <= good) return 100.0;
    if (x >= bad) return 0.0;
    return 100.0 * (1.0 - (x-good)/(bad-good));
}
static double score_high(double x, double bad, double good)
{
    if (x <= bad) return 0.0;
    if (x >= good) return 100.0;
    return 100.0 * (x-bad)/(good-bad);
}
static double score_band(double x, double bad_lo, double good_lo,
                         double good_hi, double bad_hi)
{
    if (x <= bad_lo || x >= bad_hi) return 0.0;
    if (x >= good_lo && x <= good_hi) return 100.0;
    if (x < good_lo) return 100.0 * (x-bad_lo)/(good_lo-bad_lo);
    return 100.0 * (bad_hi-x)/(bad_hi-good_hi);
}
static int sector_index(double angle, int n)
{
    double a = angle;
    if (a < 0) a += 2.0*M_PI;
    int i = (int)floor(a/(2.0*M_PI)*n);
    if (i < 0) i = 0;
    if (i >= n) i = n-1;
    return i;
}
static int best_circular_run(const int *v, int n)
{
    int best=0, cur=0;
    for (int i=0; i<2*n; i++) {
        if (v[i%n]) {
            if (++cur > best) best=cur;
            if (best >= n) return n;
        } else cur=0;
    }
    return best > n ? n : best;
}
static int circular_groups(const int *v, int n)
{
    int any=0, all=1, groups=0;
    for (int i=0;i<n;i++) { if(v[i]) any=1; else all=0; }
    if (!any) return 0;
    if (all) return 1;
    for (int i=0;i<n;i++)
        if (v[i] && !v[(i+n-1)%n]) groups++;
    return groups;
}
static double sample_height_biome(Generator *g, SurfaceNoise *sn,
                                  int bx, int bz, int *biome)
{
    int qx = bx / 4;
    int qz = bz / 4;
    float y = 0.0f;
    int id = -1;
    if (mapApproxHeight(&y, &id, g, sn, qx, qz, 1, 1) != 0)
        return NAN;
    if (biome) *biome = id;
    return (double)y;
}
static int oceanic(int id)
{
    return isOceanic(id);
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

    Stats core, core_sec[SECTORS], edge_sec[SECTORS];
    Stats trans, trans_sec[SECTORS];
    Stats ring, ring_sec[SECTORS];
    long outer_n[OUTER_SECTORS] = {0};
    long outer_land_n[OUTER_SECTORS] = {0};
    long outer_total=0, outer_land_total=0;

    stats_init(&core); stats_init(&trans); stats_init(&ring);
    for (int i=0;i<SECTORS;i++) {
        stats_init(&core_sec[i]); stats_init(&edge_sec[i]);
        stats_init(&trans_sec[i]); stats_init(&ring_sec[i]);
    }

    double max_ring = -1e30;

    for (int z=-MAX_RADIUS; z<=MAX_RADIUS; z+=STEP_BLOCKS) {
        for (int x=-MAX_RADIUS; x<=MAX_RADIUS; x+=STEP_BLOCKS) {
            double d = hypot((double)x,(double)z);
            if (d > MAX_RADIUS) continue;

            int biome=-1;
            double h = sample_height_biome(&g,&sn,x,z,&biome);
            if (!isfinite(h)) continue;
            int s8 = sector_index(atan2((double)z,(double)x),SECTORS);

            if (d <= CORE_RADIUS) {
                stats_add(&core,h);
                stats_add(&core_sec[s8],h);
            }
            if (d >= CORE_EDGE_INNER && d <= CORE_EDGE_OUTER)
                stats_add(&edge_sec[s8],h);

            if (d > TRANS_INNER && d <= TRANS_OUTER) {
                stats_add(&trans,h);
                stats_add(&trans_sec[s8],h);
            }
            if (d > RING_INNER && d <= RING_OUTER) {
                stats_add(&ring,h);
                stats_add(&ring_sec[s8],h);
                if (h > max_ring) max_ring=h;
            }
            if (d > OUTER_INNER && d <= OUTER_OUTER) {
                int s16=sector_index(atan2((double)z,(double)x),OUTER_SECTORS);
                outer_n[s16]++; outer_total++;
                if (!oceanic(biome)) {
                    outer_land_n[s16]++; outer_land_total++;
                }
            }
        }
    }

    double cm=stats_mean(&core), cs=stats_sd(&core);
    double cr=(core.n ? core.maxv-core.minv : 999.0);
    double tm=stats_mean(&trans);
    double rm=stats_mean(&ring), rs=stats_sd(&ring);

    r.core_mean=cm; r.core_sd=cs; r.core_range=cr;
    r.transition_mean=tm; r.transition_rise=tm-cm;
    r.ring_mean=rm; r.ring_sd=rs;
    r.max_ring_height=max_ring;
    r.max_ring_rise=max_ring-cm;

    /* 1) Bassin central : petit et reellement constructible.
       On combine dispersion et amplitude pour eviter qu'une seule bosse soit masquee. */
    double flat_sd = score_low(cs, 4.5, 12.0);
    double flat_range = score_low(cr, 18.0, 48.0);
    r.core_flat = 0.62*flat_sd + 0.38*flat_range;

    /* 2) Circularite du bassin.
       Les 8 directions doivent rester proches de la moyenne du centre.
       On regarde surtout l'anneau 120-175, qui correspond a la limite souhaitee. */
    int calm_edges=0;
    double worst_edge_dev=0.0;
    for (int i=0;i<SECTORS;i++) {
        if (!edge_sec[i].n) continue;
        double dev=fabs(stats_mean(&edge_sec[i])-cm);
        if (dev>worst_edge_dev) worst_edge_dev=dev;
        if (dev <= 8.0 && stats_sd(&edge_sec[i]) <= 9.0) calm_edges++;
    }
    double calm_score = score_high((double)calm_edges,4.0,7.0);
    double worst_score = score_low(worst_edge_dev,8.0,22.0);
    r.core_circular = 0.70*calm_score + 0.30*worst_score;

    /* 3) Relief qui commence vite : 150-250 doit deja monter.
       8-20 blocs de hausse moyenne est la bonne zone ; trop peu = montagnes lointaines. */
    double trans_rise=tm-cm;
    double rise_score=score_band(trans_rise,2.0,8.0,22.0,38.0);
    int rising_dirs=0;
    for (int i=0;i<SECTORS;i++) {
        if (trans_sec[i].n && stats_mean(&trans_sec[i])-cm >= 7.0)
            rising_dirs++;
    }
    double dirs_score=score_high((double)rising_dirs,3.0,6.0);
    r.close_slope=0.60*rise_score+0.40*dirs_score;

    /* 4) Couronne 200-450 : chaque secteur doit avoir un vrai relief local.
       On ne recompense plus la simple variance globale. */
    int good[SECTORS]={0};
    double sector_strength_sum=0.0;
    for (int i=0;i<SECTORS;i++) {
        if (!ring_sec[i].n) continue;
        double meanrise=stats_mean(&ring_sec[i])-cm;
        double maxrise=ring_sec[i].maxv-cm;
        double local=0.55*score_high(meanrise,7.0,20.0)
                    +0.45*score_high(maxrise,22.0,50.0);
        sector_strength_sum += local;
        if (meanrise >= 10.0 && maxrise >= 30.0) good[i]=1;
    }
    r.ring_strength=sector_strength_sum/(double)SECTORS;

    int ngood=0;
    for (int i=0;i<SECTORS;i++) ngood+=good[i];
    r.ring_good_sectors=ngood;
    r.ring_best_good_run=best_circular_run(good,SECTORS);

    /* 5) Enceinte : ideal 6-7/8. 8/8 reste bon mais pas optimal. */
    if (ngood==7) r.enclosure=100.0;
    else if (ngood==6) r.enclosure=96.0;
    else if (ngood==8) r.enclosure=88.0;
    else if (ngood==5) r.enclosure=76.0;
    else if (ngood==4) r.enclosure=50.0;
    else if (ngood==3) r.enclosure=25.0;
    else r.enclosure=0.0;

    /* 6) Forme des ouvertures : 1-2 groupes, idealement 1-2 secteurs au total.
       Penalise les reliefs disposes en deux gros massifs opposes. */
    int open[SECTORS];
    for (int i=0;i<SECTORS;i++) open[i]=!good[i];
    int groups=circular_groups(open,SECTORS);
    int openrun=best_circular_run(open,SECTORS);
    r.opening_groups=groups;
    r.largest_opening_run=openrun;

    double group_score;
    if (groups==1) group_score=100.0;
    else if (groups==2) group_score=85.0;
    else if (groups==0) group_score=72.0; /* totalement ferme */
    else if (groups==3) group_score=42.0;
    else group_score=10.0;

    double width_score;
    if (openrun==1) width_score=100.0;
    else if (openrun==2) width_score=92.0;
    else if (openrun==0) width_score=72.0;
    else if (openrun==3) width_score=55.0;
    else width_score=15.0;
    r.openings=0.55*group_score+0.45*width_score;

    /* 7) Troisieme couronne / zoo : on mesure surtout la terre disponible,
       pas sa platitude. Ocean sur un cote est acceptable, mais pas partout. */
    int land_sector[OUTER_SECTORS]={0};
    int land_sectors=0;
    for (int i=0;i<OUTER_SECTORS;i++) {
        double f = outer_n[i] ? (double)outer_land_n[i]/outer_n[i] : 0.0;
        if (f >= 0.70) { land_sector[i]=1; land_sectors++; }
    }
    int landrun=best_circular_run(land_sector,OUTER_SECTORS);
    double landfrac=outer_total ? (double)outer_land_total/outer_total : 0.0;
    r.outer_land_fraction=landfrac;
    r.outer_land_sectors=land_sectors;
    r.outer_best_land_run=landrun;

    double frac_score=score_high(landfrac,0.50,0.78);
    double run_score=score_high((double)landrun,5.0,10.0);
    double count_score=score_high((double)land_sectors,8.0,13.0);
    r.outer_land=0.50*frac_score+0.30*run_score+0.20*count_score;

    /* 8) Sommet remarquable : bonus QoL uniquement. */
    r.peak_bonus=score_band(max_ring-cm,32.0,55.0,105.0,155.0);

    /* Score structurel. Le bassin + couronne + enceinte dominent.
       Le sommet ne peut pas sauver une mauvaise geometrie. */
    r.score =
          0.20*r.core_flat
        + 0.12*r.core_circular
        + 0.15*r.close_slope
        + 0.18*r.ring_strength
        + 0.15*r.enclosure
        + 0.08*r.openings
        + 0.10*r.outer_land
        + 0.02*r.peak_bonus;

    /* Gates souples : un mauvais coeur ou une couronne trop incomplete
       doit chuter, meme si les autres metriques sont bonnes. */
    if (r.core_flat < 35.0) r.score *= 0.72;
    if (r.ring_good_sectors < 4) r.score *= 0.70;
    if (r.outer_land_fraction < 0.42) r.score *= 0.80;

    return r;
}

static int cmp_desc(const void *a,const void *b)
{
    const Result *x=(const Result*)a,*y=(const Result*)b;
    if (x->score<y->score) return 1;
    if (x->score>y->score) return -1;
    return x->seed<y->seed ? -1 : x->seed>y->seed;
}

int main(int argc,char **argv)
{
    if (argc<3) {
        fprintf(stderr,"Usage: %s seeds.txt classement_v2.csv\n",argv[0]);
        return 1;
    }
    FILE *f=fopen(argv[1],"r");
    if(!f){perror("seeds");return 1;}

    size_t cap=10000,n=0;
    uint64_t *seeds=malloc(cap*sizeof(*seeds));
    char line[256];
    while(fgets(line,sizeof(line),f)){
        char *p=line,*end;
        while(*p==' '||*p=='\t'||*p=='\r'||*p=='\n')p++;
        if(!*p||*p=='#')continue;
        unsigned long long v=strtoull(p,&end,10);
        if(end==p)continue;
        if(n==cap){cap*=2;seeds=realloc(seeds,cap*sizeof(*seeds));}
        seeds[n++]=(uint64_t)v;
    }
    fclose(f);
    if(!n){fprintf(stderr,"Aucune seed.\n");free(seeds);return 1;}

    Result *res=malloc(n*sizeof(*res));
    if(!res){free(seeds);return 1;}

    fprintf(stderr,"Scoring V2 de %zu seeds...\n",n);
    for(size_t i=0;i<n;i++){
        res[i]=score_seed(seeds[i]);
        if((i+1)%25==0||i+1==n){
            fprintf(stderr,"\r%zu / %zu",i+1,n);fflush(stderr);
        }
    }
    fprintf(stderr,"\nTri...\n");
    qsort(res,n,sizeof(*res),cmp_desc);

    FILE *o=fopen(argv[2],"w");
    if(!o){perror("csv");free(res);free(seeds);return 1;}
    fprintf(o,
      "rank,seed,score,core_flat,core_circular,close_slope,ring_strength,"
      "enclosure,openings,outer_land,peak_bonus,core_mean,core_sd,core_range,"
      "transition_mean,transition_rise,ring_mean,ring_sd,ring_good_sectors,"
      "ring_best_good_run,opening_groups,largest_opening_run,"
      "outer_land_fraction,outer_land_sectors,outer_best_land_run,"
      "max_ring_height,max_ring_rise\n");

    for(size_t i=0;i<n;i++){
        Result *q=&res[i];
        fprintf(o,
          "%zu,%llu,%.3f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,"
          "%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d,%d,%d,"
          "%.4f,%d,%d,%.3f,%.3f\n",
          i+1,(unsigned long long)q->seed,q->score,
          q->core_flat,q->core_circular,q->close_slope,q->ring_strength,
          q->enclosure,q->openings,q->outer_land,q->peak_bonus,
          q->core_mean,q->core_sd,q->core_range,
          q->transition_mean,q->transition_rise,q->ring_mean,q->ring_sd,
          q->ring_good_sectors,q->ring_best_good_run,
          q->opening_groups,q->largest_opening_run,
          q->outer_land_fraction,q->outer_land_sectors,q->outer_best_land_run,
          q->max_ring_height,q->max_ring_rise);
    }
    fclose(o);

    fprintf(stderr,"Termine: %s\nTop 25:\n",argv[2]);
    for(size_t i=0;i<n&&i<25;i++)
        fprintf(stderr,"%2zu. seed=%llu score=%.2f core=%.1f ring=%.1f enc=%.1f outer=%.1f\n",
          i+1,(unsigned long long)res[i].seed,res[i].score,
          res[i].core_flat,res[i].ring_strength,res[i].enclosure,res[i].outer_land);

    free(res);free(seeds);
    return 0;
}
