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
          R_outside_ins_or_del_seq    mh koh_mh        COSMIC_83
                             <int> <int>  <int>           <char>
       1:                       20     0      0        DEL:C:1:1
       2:                       21     0      0        DEL:T:1:0
       3:                       17     0      0        INS:T:1:4
       4:                       20     0      0        DEL:C:1:1
       5:                       20     0      0        DEL:T:1:1
       6:                       21     0      0        DEL:T:1:0
       7:                       21     0      0        INS:T:1:0
       8:                       21     0      0        INS:C:1:0
       9:                       10     0      0 INS:repeats:2:5+
      10:                       18     0      0        DEL:C:1:3
      11:                       19     0      0        DEL:T:1:2
      12:                       20     0      0        DEL:C:1:1
      13:                       17     0      0        INS:T:1:4
      14:                        0     2      2      DEL:MH:5+:2
      15:                       19     0      0        DEL:T:1:2
      16:                        0     3      3       DEL:MH:4:3
      17:                       18     0      0        DEL:T:1:3
      18:                       20     0      0        DEL:C:1:1
      19:                       18     0      0        DEL:C:1:3
      20:                        0     0      0  DEL:repeats:2:0
      21:                        1     0      0  DEL:repeats:3:1
      22:                       16     0      0       INS:T:1:5+
      23:                        0     1      1      DEL:MH:5+:1
      24:                        0     3      3      DEL:MH:5+:3
      25:                       14     0      0       INS:C:1:5+
      26:                       20     0      0        DEL:C:1:1
      27:                       20     0      0        INS:T:1:1
      28:                       18     0      0        DEL:C:1:3
      29:                       19     0      0        DEL:T:1:2
      30:                       20     0      0        INS:T:1:1
      31:                       19     0      0        DEL:T:1:2
      32:                       21     0      0        DEL:T:1:0
      33:                       17     0      0        INS:T:1:4
      34:                       20     0      0        INS:T:1:1
      35:                       15     0      0       INS:T:1:5+
      36:                       18     0      0        INS:T:1:3
      37:                       21     0      0        DEL:C:1:0
      38:                       20     0      0        INS:T:1:1
      39:                        1     0      0       INS:T:1:5+
      40:                       21     0      0        DEL:C:1:0
      41:                       17     0      0        INS:T:1:4
      42:                       21     0      0        DEL:C:1:0
      43:                       20     0      0        INS:T:1:1
      44:                        0     2      2      DEL:MH:5+:2
      45:                       18     0      0        DEL:C:1:3
      46:                       12     0      0       INS:T:1:5+
      47:                       21     0      0        DEL:T:1:0
      48:                       21     0      0        DEL:T:1:0
      49:                       20     0      0        INS:T:1:1
      50:                       17     0      0        INS:T:1:4
      51:                       20     0      0        DEL:C:1:1
      52:                       19     0      0        DEL:C:1:2
      53:                       18     0      0        DEL:T:1:3
      54:                       21     0      0        DEL:C:1:0
      55:                       14     0      0       INS:T:1:5+
      56:                        1     0      0  DEL:repeats:2:1
      57:                       20     0      0        DEL:T:1:1
      58:                       18     0      0        INS:T:1:3
      59:                        0     1      1  INS:repeats:3:0
      60:                        0     0      0  DEL:repeats:2:0
          R_outside_ins_or_del_seq    mh koh_mh        COSMIC_83
                         Koh_89        Koh_476   prev_COSMIC_83
                         <char>         <char>           <char>
       1:          [Del(C):R2]A  G[Del(C):R2]A        DEL:C:1:1
       2:     A[Del(T):R(1,4)]G  A[Del(T):R1]G        DEL:T:1:0
       3:     C[Ins(T):R(0,4)]A  C[Ins(T):R4]A        INS:T:1:4
       4:          [Del(C):R2]T  A[Del(C):R2]T        DEL:C:1:1
       5:     A[Del(T):R(1,4)]A  A[Del(T):R2]A        DEL:T:1:1
       6:     A[Del(T):R(1,4)]C  A[Del(T):R1]C        DEL:T:1:0
       7:     G[Ins(T):R(0,4)]C  G[Ins(T):R0]C        INS:T:1:0
       8:         Ins(C):R(0,3)  G[Ins(C):R0]A        INS:C:1:0
       9:         Ins(2,):R(5,) Ins2:U2:R(5,9) INS:repeats:2:5+
      10:      [Del(C):R(4,5)]T  A[Del(C):R4]T        DEL:C:1:3
      11:     A[Del(T):R(1,4)]A  A[Del(T):R3]A        DEL:T:1:2
      12:          [Del(C):R2]T  T[Del(C):R2]T        DEL:C:1:1
      13:     G[Ins(T):R(0,4)]C  G[Ins(T):R4]C        INS:T:1:4
      14:            Del(6,):M2     Del(7,):M2      DEL:MH:5+:2
      15:     A[Del(T):R(1,4)]G  A[Del(T):R3]G        DEL:T:1:2
      16:       Del(4,5):M(3,4)        Del4:M3       DEL:MH:4:3
      17:     C[Del(T):R(1,4)]C  C[Del(T):R4]C        DEL:T:1:3
      18:          [Del(C):R2]T  A[Del(C):R2]T        DEL:C:1:1
      19:      [Del(C):R(4,5)]A  A[Del(C):R4]A        DEL:C:1:3
      20: Del(2,):U(1,2):R(2,4)     Del2:U1:R2  DEL:repeats:2:0
      21:      Del(3,):U(3,):R2     Del3:U3:R2  DEL:repeats:3:1
      22:     G[Ins(T):R(5,7)]G  G[Ins(T):R5]G       INS:T:1:5+
      23:            Del(6,):M1     Del(7,):M1      DEL:MH:5+:1
      24:            Del(6,):M3        Del6:M3      DEL:MH:5+:3
      25:          Ins(C):R(7,)  A[Ins(C):R7]A       INS:C:1:5+
      26:          [Del(C):R2]A  G[Del(C):R2]A        DEL:C:1:1
      27:     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      28:      [Del(C):R(4,5)]A  A[Del(C):R4]A        DEL:C:1:3
      29:     C[Del(T):R(1,4)]G  C[Del(T):R3]G        DEL:T:1:2
      30:     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      31:     A[Del(T):R(1,4)]A  A[Del(T):R3]A        DEL:T:1:2
      32:     G[Del(T):R(1,4)]A  G[Del(T):R1]A        DEL:T:1:0
      33:     C[Ins(T):R(0,4)]G  C[Ins(T):R4]G        INS:T:1:4
      34:     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      35:     G[Ins(T):R(5,7)]G  G[Ins(T):R6]G       INS:T:1:5+
      36:     G[Ins(T):R(0,4)]C  G[Ins(T):R3]C        INS:T:1:3
      37:          [Del(C):R1]A  G[Del(C):R1]A        DEL:C:1:0
      38:     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      39:      T[Ins(T):R(8,)]A T[Ins(T):R13]A       INS:T:1:5+
      40:          [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      41:     G[Ins(T):R(0,4)]A  G[Ins(T):R4]A        INS:T:1:4
      42:          [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      43:     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      44:            Del(6,):M2     Del(7,):M2      DEL:MH:5+:2
      45:      [Del(C):R(4,5)]T  T[Del(C):R4]T        DEL:C:1:3
      46:      G[Ins(T):R(8,)]C  G[Ins(T):R9]C       INS:T:1:5+
      47:     C[Del(T):R(1,4)]G  C[Del(T):R1]G        DEL:T:1:0
      48:     C[Del(T):R(1,4)]C  C[Del(T):R1]C        DEL:T:1:0
      49:     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      50:     G[Ins(T):R(0,4)]A  G[Ins(T):R4]A        INS:T:1:4
      51:          [Del(C):R2]T  G[Del(C):R2]T        DEL:C:1:1
      52:      [Del(C):R(1,5)]G  A[Del(C):R3]G        DEL:C:1:2
      53:     A[Del(T):R(1,4)]G  A[Del(T):R4]G        DEL:T:1:3
      54:          [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      55:     A[Ins(T):R(5,7)]A  A[Ins(T):R7]A       INS:T:1:5+
      56: Del(2,):U(1,2):R(2,4)     Del2:U2:R2  DEL:repeats:2:1
      57:     A[Del(T):R(1,4)]G  A[Del(T):R2]G        DEL:T:1:1
      58:     A[Ins(T):R(0,4)]G  A[Ins(T):R3]G        INS:T:1:3
      59:           Ins(2,4):R0     Ins(2,4):M  INS:repeats:3:0
      60: Del(2,):U(1,2):R(2,4)     Del2:U1:R2  DEL:repeats:2:0
                         Koh_89        Koh_476   prev_COSMIC_83

