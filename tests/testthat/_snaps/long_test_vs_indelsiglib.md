# test_long_test_vs_indelsiglib

    Code
      xx[zz, ]
    Output
         CHROM      POS REF ALT                                 seq.context
      82  chr3 89101431   A  AC CACCCCACAGAAGTGAGATGAACCCCCCCCCCAACTCAAAATT
         seq.context.width Koh89.annotate.class Koh476.annotate.class ins_or_del pre
      82                21                             A[Ins(C):R10]A          i   A
         ins_or_del_seq post L U_seq U U_seq_count_in_indel_seq
      82              C    A 1     C 1                        1
         indel_str_count_in_ref  R R_outside_ins_or_del_seq mh unit unit_length
      82                     10 10                       11  0    C           1
         internal_rep internal_reps spacer spacer_length prime3_rep prime3_reps
      82                          0                    0 CCCCCCCCCC          10
         original_reps  COSMIC_83       Koh_89        Koh_476 prev_COSMIC_83
      82            10 INS:C:1:5+ Ins(C):R(7,) A[Ins(C):R10]A     INS:C:1:5+
         koh_orig_edited
      82                

