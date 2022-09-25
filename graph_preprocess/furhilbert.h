//
// Created by atrostan on 14/09/22.
//

#ifndef GRAPH_PREPROCESS_FURHILBERT_H
#define GRAPH_PREPROCESS_FURHILBERT_H

/*
 * File:   furhilbert.h
 * Author: Christian BÃ¶hm, Martin Perdacher, and Claudia Plant
 *
 * Description: a very efficient generation of coordinate pairs following a new
 * variant of the Hilbert curve, called FUR-Hilbert (Fast UnRestricted). The
 * method is described in the following literature which must be cited in any
 * scientific paper that uses FUR-Hilbert:
 *
@inproceedings{furhilbert,
  author    = {Christian B{\"{o}}hm and
               Martin Perdacher and
               Claudia Plant},
  title     = {Cache-oblivious Loops Based on a Novel Space-filling Curve},
  booktitle = {IEEE Int. Conference on Big Data},
  year      = {2016},
}
 *
 * Defined is a pair of preprocessor macros (FUR_HILBERT_FOR, FUR_HILBERT_END)
 * between which the intended loop body is enclosed. Two iterator variables
 * (the name of which can be freely selected) are updated in each loop iteration
 * following a Hilbert-like style, i.e. exactly one of the variables is incremented
 * or decremented by 1 in each iteration while generally following a fair,
 * bisecting scheme. In contrast to the original Hilbert curve which is restricted
 * to squares with a side length of a power of 2, FUR-Hilbert fills any arbitrary
 * rectangle. The overhead per loop iteration is constant (in contrast to O(log n)
 * for previous methods). The purpose is cache-obliviousness: objects indexed by
 * the iterator variables have a high chance to remain in LRU caches of arbitrary
 * size.
 *
 * Usage: as demonstrated in the following example:
 *
int i, j;
FUR_HILBERT_FOR(i, j, 3, 6, 2, 10){
    printf("%d %d\n", i, j) ;
} FUR_HILBERT_END(i, j);
 *
 * This will output all pairs $(i,j) \in \{3,4,5\} \times \{2, ..., 9\}$.
 * As usual in C-like languages, the upper bounds (6 and 10) are excluded.
 * Changes of upper and lower bounds in the loop body are not allowed (ignored).
 */

#ifndef FURHILBERT_H
#define	FURHILBERT_H

extern unsigned long HILLOOP_nanoprog[9][9][4][2];

#define FUR_HILBERT_FOR(i,j,imin,imax,jmin,jmax) {\
    int HILLOOP_imin = (imin);\
    int HILLOOP_imax = (imax);\
    int HILLOOP_jmin = (jmin);\
    int HILLOOP_jmax = (jmax);\
    (i) = HILLOOP_imin;\
    (j) = HILLOOP_jmin;\
    int HILLOOP_idiff = HILLOOP_imax - HILLOOP_imin;\
    int HILLOOP_jdiff = HILLOOP_jmax - HILLOOP_jmin;\
    if (HILLOOP_idiff > 0 && HILLOOP_jdiff > 0) {\
        int HILLOOP_iceil = HILLOOP_idiff - 1;\
        HILLOOP_iceil |= HILLOOP_iceil >> 1;\
        HILLOOP_iceil |= HILLOOP_iceil >> 2;\
        HILLOOP_iceil |= HILLOOP_iceil >> 4;\
        HILLOOP_iceil |= HILLOOP_iceil >> 8;\
        HILLOOP_iceil |= HILLOOP_iceil >> 16;\
        HILLOOP_iceil++;\
        if (HILLOOP_iceil < 8)\
            HILLOOP_iceil = 8;\
        int HILLOOP_jceil = HILLOOP_jdiff - 1;\
        HILLOOP_jceil |= HILLOOP_jceil >> 1;\
        HILLOOP_jceil |= HILLOOP_jceil >> 2;\
        HILLOOP_jceil |= HILLOOP_jceil >> 4;\
        HILLOOP_jceil |= HILLOOP_jceil >> 8;\
        HILLOOP_jceil |= HILLOOP_jceil >> 16;\
        HILLOOP_jceil++;\
        if (HILLOOP_jceil < 8)\
            HILLOOP_jceil = 8;\
        while ((i) < HILLOOP_imax && (j) < HILLOOP_jmax) {\
            int HILLOOP_icur, HILLOOP_jcur, HILLOOP_c, HILLOOP_base, HILLOOP_icase, HILLOOP_jcase;\
            unsigned long long HILLOOP_stop, HILLOOP_hilbert = 0ull;\
            if (HILLOOP_idiff > HILLOOP_jdiff) {\
                HILLOOP_jcur = HILLOOP_jdiff;\
                HILLOOP_icur = HILLOOP_imax - (i);\
                HILLOOP_base = HILLOOP_jceil / 8;\
                HILLOOP_stop = HILLOOP_base * HILLOOP_base;\
                if (HILLOOP_icur >= 2 * HILLOOP_jceil) {\
                    HILLOOP_icur = HILLOOP_jceil;\
                    HILLOOP_c = 3;\
                    HILLOOP_icase = 0;\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                } else if (HILLOOP_icur > HILLOOP_jceil) {\
                    HILLOOP_icur = (HILLOOP_icur + 3) / 4 * 2;\
                    HILLOOP_c = 3;\
                    HILLOOP_icase = 0;\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                } else {\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                    HILLOOP_c = 3 - (HILLOOP_icase & ~HILLOOP_jcase);\
                    HILLOOP_icase += HILLOOP_icase & HILLOOP_jcase;\
                }\
            } else {\
                HILLOOP_icur = HILLOOP_idiff;\
                HILLOOP_jcur = HILLOOP_jmax - (j);\
                HILLOOP_base = HILLOOP_iceil / 8;\
                HILLOOP_stop = HILLOOP_base * HILLOOP_base;\
                if (HILLOOP_jcur >= 2 * HILLOOP_iceil) {\
                    HILLOOP_jcur = HILLOOP_iceil;\
                    HILLOOP_c = 2;\
                    HILLOOP_jcase = 0;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                } else if (HILLOOP_jcur > HILLOOP_iceil) {\
                    HILLOOP_jcur = (HILLOOP_jcur + 3) / 4 * 2;\
                    HILLOOP_c = 2;\
                    HILLOOP_jcase = 0;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                } else {\
                    HILLOOP_jcase = HILLOOP_jcur & 1;\
                    HILLOOP_icase = HILLOOP_icur & 1;\
                    HILLOOP_c = 2 + (HILLOOP_jcase & ~HILLOOP_icase);\
                    HILLOOP_jcase += HILLOOP_icase & HILLOOP_jcase;\
                }\
            }\
            int HILLOOP_i57 = 5 + 2 * (HILLOOP_icur > HILLOOP_base * 6);\
            if(HILLOOP_icur<6) HILLOOP_i57 = HILLOOP_icur ;\
            int HILLOOP_j57 = 5 + 2 * (HILLOOP_jcur > HILLOOP_base * 6);\
            if(HILLOOP_jcur<6) HILLOOP_j57 = HILLOOP_jcur ;\
            int HILLOOP_isize, HILLOOP_jsize;\
            HILLOOP_c ^= ((HILLOOP_base & 0x55555555) == 0) ;\
            int HILLOOP_I = 0;\
            int HILLOOP_J = 0;\
            while (HILLOOP_hilbert < HILLOOP_stop) {\
                if (HILLOOP_icase)\
                    if (HILLOOP_icase == 2) {\
                        int HILLOOP_v = HILLOOP_base - 1 - HILLOOP_J;\
                        HILLOOP_v |= HILLOOP_v >> 1;\
                        HILLOOP_v |= HILLOOP_v >> 2;\
                        HILLOOP_v |= HILLOOP_v >> 4;\
                        HILLOOP_v |= HILLOOP_v >> 8;\
                        HILLOOP_v |= HILLOOP_v >> 16;\
                        HILLOOP_v = (HILLOOP_v + 1) / 2;\
                        HILLOOP_isize = (HILLOOP_I == HILLOOP_v ? HILLOOP_i57 : (HILLOOP_I < HILLOOP_v ?\
                                (HILLOOP_I + 1) * (HILLOOP_icur-HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                                - HILLOOP_I * (HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2 :\
                                HILLOOP_I * (HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                                - (HILLOOP_I - 1) * (HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2));\
                    } else HILLOOP_isize = HILLOOP_I == HILLOOP_base - 1 ? HILLOOP_i57 :\
                            (HILLOOP_I + 1) * (HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2\
                        - HILLOOP_I * (HILLOOP_icur - HILLOOP_i57) / (HILLOOP_base - 1) / 2 * 2;\
                else HILLOOP_isize = (HILLOOP_I + 1) * HILLOOP_icur / HILLOOP_base / 2 * 2\
                        - HILLOOP_I * HILLOOP_icur / HILLOOP_base / 2 * 2;\
                if (HILLOOP_jcase)\
                    if (HILLOOP_jcase == 2) {\
                        int HILLOOP_v = HILLOOP_base - 1 - HILLOOP_I;\
                        HILLOOP_v |= HILLOOP_v >> 1;\
                        HILLOOP_v |= HILLOOP_v >> 2;\
                        HILLOOP_v |= HILLOOP_v >> 4;\
                        HILLOOP_v |= HILLOOP_v >> 8;\
                        HILLOOP_v |= HILLOOP_v >> 16;\
                        HILLOOP_v = (HILLOOP_v + 1) / 2;\
                        HILLOOP_jsize = (HILLOOP_J == HILLOOP_v ? HILLOOP_j57 : (HILLOOP_J < HILLOOP_v ?\
                                (HILLOOP_J + 1) * (HILLOOP_jcur-HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                                - HILLOOP_J * (HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2 :\
                                HILLOOP_J * (HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                                - (HILLOOP_J - 1) * (HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2));\
                    } else HILLOOP_jsize = HILLOOP_J == HILLOOP_base - 1 ? HILLOOP_j57 :\
                            (HILLOOP_J + 1) * (HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2\
                        - HILLOOP_J * (HILLOOP_jcur - HILLOOP_j57) / (HILLOOP_base - 1) / 2 * 2;\
                else HILLOOP_jsize = (HILLOOP_J + 1) * HILLOOP_jcur / HILLOOP_base / 2 * 2\
                        - HILLOOP_J * HILLOOP_jcur / HILLOOP_base / 2 * 2;\
                unsigned long long HILLOOP_nanoH = HILLOOP_nanoprog[HILLOOP_isize][HILLOOP_jsize][HILLOOP_c][0] ;\
                unsigned long long HILLOOP_nanoL = HILLOOP_nanoprog[HILLOOP_isize][HILLOOP_jsize][HILLOOP_c][1] ;\
                for(;;){


# define FUR_HILBERT_END(i,j)\
                    if(HILLOOP_nanoL == 1)\
                        break;\
                    int HILLOOP_d = HILLOOP_nanoL & 1 | HILLOOP_nanoH & 2 ;\
                    (i) += (HILLOOP_d - 2) % 2 ;\
                    (j) += (HILLOOP_d - 1) % 2 ;\
                    HILLOOP_nanoL /= 2 ;\
                    HILLOOP_nanoH /= 2 ;\
                }\
                HILLOOP_hilbert++;\
                unsigned long long HILLOOP_l = HILLOOP_hilbert & -HILLOOP_hilbert;\
                HILLOOP_l = (HILLOOP_l + HILLOOP_l/2) & 0x5555555555555555ull;\
                int HILLOOP_a = (HILLOOP_hilbert / HILLOOP_l) & 3 ;\
                int HILLOOP_isOdd = (HILLOOP_l & 0xCCCCCCCCCCCCCCCCull) != 0;\
                HILLOOP_c ^= 3 * ((HILLOOP_a==3) != HILLOOP_isOdd);\
                HILLOOP_I += (HILLOOP_c - 2) % 2 ;\
                (i) += (HILLOOP_c - 2) % 2 ;\
                HILLOOP_J += (HILLOOP_c - 1) % 2 ;\
                (j) += (HILLOOP_c - 1) % 2 ;\
                HILLOOP_c ^= ((HILLOOP_a==1) != HILLOOP_isOdd);\
            }\
        }\
    }\
}


#endif	/* FURHILBERT_H */

#endif //GRAPH_PREPROCESS_FURHILBERT_H
