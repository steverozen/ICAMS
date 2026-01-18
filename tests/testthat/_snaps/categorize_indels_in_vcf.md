# categorize_indels_in_vcf

    Code
      yy
    Output
          ins_or_del    pre ins_or_del_seq   post     L        U_seq     U
              <char> <char>         <char> <char> <int>       <char> <int>
       1:          d      G              C      A     1            C     1
       2:          d      A              T      G     1            T     1
       3:          i      C              T      A     1            T     1
       4:          d      A              C      T     1            C     1
       5:          d      A              T      A     1            T     1
       6:          d      A              T      C     1            T     1
       7:          i      G              T      C     1            T     1
       8:          i      G              C      A     1            C     1
       9:          i      G             TA      T     2           TA     2
      10:          d      A              C      T     1            C     1
      11:          d      A              T      A     1            T     1
      12:          d      T              C      T     1            C     1
      13:          i      G              T      C     1            T     1
      14:          d      A      TAGTTATAC      T     9    TAGTTATAC     9
      15:          d      A              T      G     1            T     1
      16:          d      G           CCCA      C     4         CCCA     4
      17:          d      C              T      C     1            T     1
      18:          d      A              C      T     1            C     1
      19:          d      A              C      A     1            C     1
      20:          d      T             AA      G     2            A     1
      21:          d      A            CTT      T     3          CTT     3
      22:          i      G              T      G     1            T     1
      23:          d      G   CCTGAGCTGTAC      C    12 CCTGAGCTGTAC    12
      24:          d      G         CTGTAC      C     6       CTGTAC     6
      25:          i      A              C      A     1            C     1
      26:          d      G              C      A     1            C     1
      27:          i      G              T      G     1            T     1
      28:          d      A              C      A     1            C     1
      29:          d      C              T      G     1            T     1
      30:          i      G              T      G     1            T     1
      31:          d      A              T      A     1            T     1
      32:          d      G              T      A     1            T     1
      33:          i      C              T      G     1            T     1
      34:          i      G              T      A     1            T     1
      35:          i      G              T      G     1            T     1
      36:          i      G              T      C     1            T     1
      37:          d      G              C      A     1            C     1
      38:          i      G              T      G     1            T     1
      39:          i      T              T      A     1            T     1
      40:          d      T              C      T     1            C     1
      41:          i      G              T      A     1            T     1
      42:          d      T              C      T     1            C     1
      43:          i      G              T      A     1            T     1
      44:          d      A        TGATTCT      T     7      TGATTCT     7
      45:          d      T              C      T     1            C     1
      46:          i      G              T      C     1            T     1
      47:          d      C              T      G     1            T     1
      48:          d      C              T      C     1            T     1
      49:          i      G              T      A     1            T     1
      50:          i      G              T      A     1            T     1
      51:          d      G              C      T     1            C     1
      52:          d      A              C      G     1            C     1
      53:          d      A              T      G     1            T     1
      54:          d      T              C      T     1            C     1
      55:          i      A              T      A     1            T     1
      56:          d      C             TG      T     2           TG     2
      57:          d      A              T      G     1            T     1
      58:          i      A              T      G     1            T     1
      59:          i      C            TCT      T     3          TCT     3
      60:          d      A             GG      T     2            G     1
          ins_or_del    pre ins_or_del_seq   post     L        U_seq     U
              <char> <char>         <char> <char> <int>       <char> <int>
          U_seq_count_in_indel_seq indel_str_count_in_ref     R
                             <int>                  <int> <int>
       1:                        1                      2     2
       2:                        1                      1     1
       3:                        1                      4     4
       4:                        1                      2     2
       5:                        1                      2     2
       6:                        1                      1     1
       7:                        1                      0     0
       8:                        1                      0     0
       9:                        1                     10    10
      10:                        1                      4     4
      11:                        1                      3     3
      12:                        1                      2     2
      13:                        1                      4     4
      14:                        1                      1     1
      15:                        1                      3     3
      16:                        1                      1     1
      17:                        1                      4     4
      18:                        1                      2     2
      19:                        1                      4     4
      20:                        2                      1     2
      21:                        1                      2     2
      22:                        1                      5     5
      23:                        1                      1     1
      24:                        1                      1     1
      25:                        1                      7     7
      26:                        1                      2     2
      27:                        1                      1     1
      28:                        1                      4     4
      29:                        1                      3     3
      30:                        1                      1     1
      31:                        1                      3     3
      32:                        1                      1     1
      33:                        1                      4     4
      34:                        1                      1     1
      35:                        1                      6     6
      36:                        1                      3     3
      37:                        1                      1     1
      38:                        1                      1     1
      39:                        1                     20    20
      40:                        1                      1     1
      41:                        1                      4     4
      42:                        1                      1     1
      43:                        1                      1     1
      44:                        1                      1     1
      45:                        1                      4     4
      46:                        1                      9     9
      47:                        1                      1     1
      48:                        1                      1     1
      49:                        1                      1     1
      50:                        1                      4     4
      51:                        1                      2     2
      52:                        1                      3     3
      53:                        1                      4     4
      54:                        1                      1     1
      55:                        1                      7     7
      56:                        1                      2     2
      57:                        1                      2     2
      58:                        1                      3     3
      59:                        1                      0     0
      60:                        2                      1     2
          U_seq_count_in_indel_seq indel_str_count_in_ref     R
                             <int>                  <int> <int>
          R_outside_ins_or_del_seq    mh koh_mh   unit unit_length internal_rep
                             <int> <int>  <int> <char>       <int>       <char>
       1:                       20     0      0      G           1             
       2:                       21     0      0      T           1             
       3:                       17     0      0      A           1             
       4:                       20     0      0      C           1             
       5:                       20     0      0      A           1             
       6:                       21     0      0      A           1             
       7:                       21     0      0      T           1             
       8:                       21     0      0      G           1             
       9:                       10     0      0     TA           2             
      10:                       18     0      0      C           1             
      11:                       19     0      0      A           1             
      12:                       20     0      0      G           1             
      13:                       17     0      0      T           1             
      14:                        0     2      2     TA           2             
      15:                       19     0      0      A           1             
      16:                        0     3      3      C           1           CC
      17:                       18     0      0      A           1             
      18:                       20     0      0      G           1             
      19:                       18     0      0      G           1             
      20:                        0     0      0      A           1            A
      21:                        1     0      0    CTT           3             
      22:                       16     0      0      T           1             
      23:                        0     1      1      C           1            C
      24:                        0     3      3    CTG           3             
      25:                       14     0      0      G           1             
      26:                       20     0      0      G           1             
      27:                       20     0      0      T           1             
      28:                       18     0      0      C           1             
      29:                       19     0      0      T           1             
      30:                       20     0      0      T           1             
      31:                       19     0      0      A           1             
      32:                       21     0      0      T           1             
      33:                       17     0      0      T           1             
      34:                       20     0      0      T           1             
      35:                       15     0      0      T           1             
      36:                       18     0      0      T           1             
      37:                       21     0      0      C           1             
      38:                       20     0      0      T           1             
      39:                        1     0      0      A           1             
      40:                       21     0      0      G           1             
      41:                       17     0      0      T           1             
      42:                       21     0      0      C           1             
      43:                       20     0      0      A           1             
      44:                        0     2      2     TG           2             
      45:                       18     0      0      G           1             
      46:                       12     0      0      T           1             
      47:                       21     0      0      A           1             
      48:                       21     0      0      A           1             
      49:                       20     0      0      A           1             
      50:                       17     0      0      T           1             
      51:                       20     0      0      C           1             
      52:                       19     0      0      C           1             
      53:                       18     0      0      T           1             
      54:                       21     0      0      C           1             
      55:                       14     0      0      A           1             
      56:                        1     0      0     TG           2             
      57:                       20     0      0      A           1             
      58:                       18     0      0      A           1             
      59:                        0     1      1      T           1             
      60:                        0     0      0      G           1            G
          R_outside_ins_or_del_seq    mh koh_mh   unit unit_length internal_rep
                             <int> <int>  <int> <char>       <int>       <char>
          internal_reps     spacer spacer_length            prime3_rep prime3_reps
                  <int>     <char>         <int>                <char>       <int>
       1:             0                        0                     G           1
       2:             0                        0                                 0
       3:             0                        0                  AAAA           4
       4:             0                        0                     C           1
       5:             0                        0                     A           1
       6:             0                        0                                 0
       7:             0                        0                                 0
       8:             0                        0                                 0
       9:             0                        0  TATATATATATATATATATA          10
      10:             0                        0                   CCC           3
      11:             0                        0                    AA           2
      12:             0                        0                     G           1
      13:             0                        0                  TTTT           4
      14:             0    GTTATAC             7                    TA           1
      15:             0                        0                    AA           2
      16:             2          A             1                  CCCC           4
      17:             0                        0                   AAA           3
      18:             0                        0                     G           1
      19:             0                        0                   GGG           3
      20:             1                        0                                 0
      21:             0                        0                   CTT           1
      22:             0                        0                 TTTTT           5
      23:             1 TGAGCTGTAC            10                     C           1
      24:             0        TAC             3                   CTG           1
      25:             0                        0               GGGGGGG           7
      26:             0                        0                     G           1
      27:             0                        0                     T           1
      28:             0                        0                   CCC           3
      29:             0                        0                    TT           2
      30:             0                        0                     T           1
      31:             0                        0                    AA           2
      32:             0                        0                                 0
      33:             0                        0                  TTTT           4
      34:             0                        0                     T           1
      35:             0                        0                TTTTTT           6
      36:             0                        0                   TTT           3
      37:             0                        0                                 0
      38:             0                        0                     T           1
      39:             0                        0 AAAAAAAAAAAAAAAAAAAAA          21
      40:             0                        0                                 0
      41:             0                        0                  TTTT           4
      42:             0                        0                                 0
      43:             0                        0                     A           1
      44:             0      ATTCT             5                    TG           1
      45:             0                        0                   GGG           3
      46:             0                        0             TTTTTTTTT           9
      47:             0                        0                                 0
      48:             0                        0                                 0
      49:             0                        0                     A           1
      50:             0                        0                  TTTT           4
      51:             0                        0                     C           1
      52:             0                        0                    CC           2
      53:             0                        0                   TTT           3
      54:             0                        0                                 0
      55:             0                        0               AAAAAAA           7
      56:             0                        0                    TG           1
      57:             0                        0                     A           1
      58:             0                        0                   AAA           3
      59:             0         CT             2           TTTTTTTTTTT          11
      60:             1                        0                                 0
          internal_reps     spacer spacer_length            prime3_rep prime3_reps
                  <int>     <char>         <int>                <char>       <int>
          original_reps        COSMIC_83                 Koh_89        Koh_476
                  <int>           <char>                 <char>         <char>
       1:             2        DEL:C:1:1           [Del(C):R2]A  G[Del(C):R2]A
       2:             1        DEL:T:1:0      A[Del(T):R(1,4)]G  A[Del(T):R1]G
       3:             4        INS:T:1:4      C[Ins(T):R(0,4)]A  C[Ins(T):R4]A
       4:             2        DEL:C:1:1           [Del(C):R2]T  A[Del(C):R2]T
       5:             2        DEL:T:1:1      A[Del(T):R(1,4)]A  A[Del(T):R2]A
       6:             1        DEL:T:1:0      A[Del(T):R(1,4)]C  A[Del(T):R1]C
       7:             0        INS:T:1:0      G[Ins(T):R(0,4)]C  G[Ins(T):R0]C
       8:             0        INS:C:1:0          Ins(C):R(0,3)  G[Ins(C):R0]A
       9:            10 INS:repeats:2:5+          Ins(2,):R(5,) Ins2:U2:R(5,9)
      10:             4        DEL:C:1:3       [Del(C):R(4,5)]T  A[Del(C):R4]T
      11:             3        DEL:T:1:2      A[Del(T):R(1,4)]A  A[Del(T):R3]A
      12:             2        DEL:C:1:1           [Del(C):R2]T  T[Del(C):R2]T
      13:             4        INS:T:1:4      G[Ins(T):R(0,4)]C  G[Ins(T):R4]C
      14:             1      DEL:MH:5+:2             Del(6,):M2     Del(7,):M2
      15:             3        DEL:T:1:2      A[Del(T):R(1,4)]G  A[Del(T):R3]G
      16:             4       DEL:MH:4:3        Del(4,5):M(3,4)        Del4:M3
      17:             4        DEL:T:1:3      C[Del(T):R(1,4)]C  C[Del(T):R4]C
      18:             2        DEL:C:1:1           [Del(C):R2]T  A[Del(C):R2]T
      19:             4        DEL:C:1:3       [Del(C):R(4,5)]A  A[Del(C):R4]A
      20:             2  DEL:repeats:2:0            Del(2,4):R1     Del2:U1:R1
      21:             2  DEL:repeats:3:1       Del(3,):U(3,):R2     Del3:U3:R2
      22:             5       INS:T:1:5+      G[Ins(T):R(5,7)]G  G[Ins(T):R5]G
      23:             1      DEL:MH:5+:1             Del(6,):M1     Del(7,):M1
      24:             1      DEL:MH:5+:3             Del(6,):M3        Del6:M3
      25:             7       INS:C:1:5+           Ins(C):R(7,)  A[Ins(C):R7]A
      26:             2        DEL:C:1:1           [Del(C):R2]A  G[Del(C):R2]A
      27:             1        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      28:             4        DEL:C:1:3       [Del(C):R(4,5)]A  A[Del(C):R4]A
      29:             3        DEL:T:1:2      C[Del(T):R(1,4)]G  C[Del(T):R3]G
      30:             1        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      31:             3        DEL:T:1:2      A[Del(T):R(1,4)]A  A[Del(T):R3]A
      32:             1        DEL:T:1:0      G[Del(T):R(1,4)]A  G[Del(T):R1]A
      33:             4        INS:T:1:4      C[Ins(T):R(0,4)]G  C[Ins(T):R4]G
      34:             1        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      35:             6       INS:T:1:5+      G[Ins(T):R(5,7)]G  G[Ins(T):R6]G
      36:             3        INS:T:1:3      G[Ins(T):R(0,4)]C  G[Ins(T):R3]C
      37:             1        DEL:C:1:0           [Del(C):R1]A  G[Del(C):R1]A
      38:             1        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      39:            21       INS:T:1:5+       T[Ins(T):R(8,)]A T[Ins(T):R20]A
      40:             1        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      41:             4        INS:T:1:4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      42:             1        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      43:             1        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      44:             1      DEL:MH:5+:2             Del(6,):M2     Del(7,):M2
      45:             4        DEL:C:1:3       [Del(C):R(4,5)]T  T[Del(C):R4]T
      46:             9       INS:T:1:5+       G[Ins(T):R(8,)]C  G[Ins(T):R9]C
      47:             1        DEL:T:1:0      C[Del(T):R(1,4)]G  C[Del(T):R1]G
      48:             1        DEL:T:1:0      C[Del(T):R(1,4)]C  C[Del(T):R1]C
      49:             1        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      50:             4        INS:T:1:4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      51:             2        DEL:C:1:1           [Del(C):R2]T  G[Del(C):R2]T
      52:             3        DEL:C:1:2       [Del(C):R(1,5)]G  A[Del(C):R3]G
      53:             4        DEL:T:1:3      A[Del(T):R(1,4)]G  A[Del(T):R4]G
      54:             1        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      55:             7       INS:T:1:5+      A[Ins(T):R(5,7)]A  A[Ins(T):R7]A
      56:             2  DEL:repeats:2:1 Del(2,8):U(1,2):R(2,4)     Del2:U2:R2
      57:             2        DEL:T:1:1      A[Del(T):R(1,4)]G  A[Del(T):R2]G
      58:             3        INS:T:1:3      A[Ins(T):R(0,4)]G  A[Ins(T):R3]G
      59:            11  INS:repeats:3:0            Ins(2,4):R0     Ins(2,4):M
      60:             2  DEL:repeats:2:0            Del(2,4):R1     Del2:U1:R1
          original_reps        COSMIC_83                 Koh_89        Koh_476
                  <int>           <char>                 <char>         <char>
            prev_COSMIC_83
                    <char>
       1:        DEL:C:1:1
       2:        DEL:T:1:0
       3:        INS:T:1:4
       4:        DEL:C:1:1
       5:        DEL:T:1:1
       6:        DEL:T:1:0
       7:        INS:T:1:0
       8:        INS:C:1:0
       9: INS:repeats:2:5+
      10:        DEL:C:1:3
      11:        DEL:T:1:2
      12:        DEL:C:1:1
      13:        INS:T:1:4
      14:      DEL:MH:5+:2
      15:        DEL:T:1:2
      16:       DEL:MH:4:3
      17:        DEL:T:1:3
      18:        DEL:C:1:1
      19:        DEL:C:1:3
      20:  DEL:repeats:2:0
      21:  DEL:repeats:3:1
      22:       INS:T:1:5+
      23:      DEL:MH:5+:1
      24:      DEL:MH:5+:3
      25:       INS:C:1:5+
      26:        DEL:C:1:1
      27:        INS:T:1:1
      28:        DEL:C:1:3
      29:        DEL:T:1:2
      30:        INS:T:1:1
      31:        DEL:T:1:2
      32:        DEL:T:1:0
      33:        INS:T:1:4
      34:        INS:T:1:1
      35:       INS:T:1:5+
      36:        INS:T:1:3
      37:        DEL:C:1:0
      38:        INS:T:1:1
      39:       INS:T:1:5+
      40:        DEL:C:1:0
      41:        INS:T:1:4
      42:        DEL:C:1:0
      43:        INS:T:1:1
      44:      DEL:MH:5+:2
      45:        DEL:C:1:3
      46:       INS:T:1:5+
      47:        DEL:T:1:0
      48:        DEL:T:1:0
      49:        INS:T:1:1
      50:        INS:T:1:4
      51:        DEL:C:1:1
      52:        DEL:C:1:2
      53:        DEL:T:1:3
      54:        DEL:C:1:0
      55:       INS:T:1:5+
      56:  DEL:repeats:2:1
      57:        DEL:T:1:1
      58:        INS:T:1:3
      59:  INS:repeats:3:0
      60:  DEL:repeats:2:0
            prev_COSMIC_83
                    <char>

