# categorize_many_indels

    Code
      yy
    Output
          ins_or_del    pre ins_or_del_seq indel_str_count_in_ref   post    mh     R
              <char> <char>         <char>                  <int> <char> <int> <int>
       1:          d      G              C                      2      A     0     2
       2:          d      A              T                      1      G     0     1
       3:          i      C              T                      4      A     0     4
       4:          d      A              C                      2      T     0     2
       5:          d      A              T                      2      A     0     2
       6:          d      A              T                      1      C     0     1
       7:          i      G              T                      0      C     0     0
       8:          i      G              C                      0      A     0     0
       9:          i      G             TA                     10      T     0    10
      10:          d      A              C                      4      T     0     4
      11:          d      A              T                      3      A     0     3
      12:          d      T              C                      2      T     0     2
      13:          i      G              T                      4      C     0     4
      14:          d      A      TAGTTATAC                      1      T     2     1
      15:          d      A              T                      3      G     0     3
      16:          d      G           CCCA                      1      C     3     1
      17:          d      C              T                      4      C     0     4
      18:          d      A              C                      2      T     0     2
      19:          d      A              C                      4      A     0     4
      20:          d      T             AA                      1      G     0     2
      21:          d      A            CTT                      2      T     0     2
      22:          i      G              T                      5      G     0     5
      23:          d      G   CCTGAGCTGTAC                      1      C     1     1
      24:          d      G         CTGTAC                      1      C     3     1
      25:          i      A              C                      7      A     0     7
      26:          d      G              C                      2      A     0     2
      27:          i      G              T                      1      G     0     1
      28:          d      A              C                      4      A     0     4
      29:          d      C              T                      3      G     0     3
      30:          i      G              T                      1      G     0     1
      31:          d      A              T                      3      A     0     3
      32:          d      G              T                      1      A     0     1
      33:          i      C              T                      4      G     0     4
      34:          i      G              T                      1      A     0     1
      35:          i      G              T                      6      G     0     6
      36:          i      G              T                      3      C     0     3
      37:          d      G              C                      1      A     0     1
      38:          i      G              T                      1      G     0     1
      39:          i      T              T                     20      A     0    20
      40:          d      T              C                      1      T     0     1
      41:          i      G              T                      4      A     0     4
      42:          d      T              C                      1      T     0     1
      43:          i      G              T                      1      A     0     1
      44:          d      A        TGATTCT                      1      T     2     1
      45:          d      T              C                      4      T     0     4
      46:          i      G              T                      9      C     0     9
      47:          d      C              T                      1      G     0     1
      48:          d      C              T                      1      C     0     1
      49:          i      G              T                      1      A     0     1
      50:          i      G              T                      4      A     0     4
      51:          d      G              C                      2      T     0     2
      52:          d      A              C                      3      G     0     3
      53:          d      A              T                      4      G     0     4
      54:          d      T              C                      1      T     0     1
      55:          i      A              T                      7      A     0     7
      56:          d      C             TG                      2      T     0     2
      57:          d      A              T                      2      G     0     2
      58:          i      A              T                      3      G     0     3
      59:          i      C            TCT                      0      T     1     0
      60:          d      A             GG                      1      T     0     2
          ins_or_del    pre ins_or_del_seq indel_str_count_in_ref   post    mh     R
              U koh_mh        COSMIC_83                Koh_89        Koh_476
          <int>  <int>           <char>                <char>         <char>
       1:     1      0        DEL:C:1:1         G[Del(C):R2]A  G[Del(C):R2]A
       2:     1      0        DEL:T:1:0      A[Del(T):R(8,)]G  A[Del(T):R1]G
       3:     1      0        INS:T:1:4     C[Ins(T):R(0,4)]A  C[Ins(T):R4]A
       4:     1      0        DEL:C:1:1         A[Del(C):R2]T  A[Del(C):R2]T
       5:     1      0        DEL:T:1:1      A[Del(T):R(8,)]A  A[Del(T):R2]A
       6:     1      0        DEL:T:1:0      A[Del(T):R(8,)]C  A[Del(T):R1]C
       7:     1      0        INS:T:1:0     G[Ins(T):R(0,4)]C  G[Ins(T):R0]C
       8:     1      0        INS:C:1:0         Ins(C):R(0,3)  G[Ins(C):R0]A
       9:     2      0 INS:repeats:2:5+         Ins(2,):R(5,)  Ins2:U2:R(5,)
      10:     1      0        DEL:C:1:3     A[Del(C):R(4,5)]T  A[Del(C):R4]T
      11:     1      0        DEL:T:1:2      A[Del(T):R(8,)]A  A[Del(T):R3]A
      12:     1      0        DEL:C:1:1         T[Del(C):R2]T  T[Del(C):R2]T
      13:     1      0        INS:T:1:4     G[Ins(T):R(0,4)]C  G[Ins(T):R4]C
      14:     9      2      DEL:MH:5+:2            Del(6,):M2     Del(7,):M2
      15:     1      0        DEL:T:1:2      A[Del(T):R(8,)]G  A[Del(T):R3]G
      16:     4      3       DEL:MH:4:3       del(2,3):M(3,4)        Del4:M3
      17:     1      0        DEL:T:1:3      C[Del(T):R(8,)]C  C[Del(T):R4]C
      18:     1      0        DEL:C:1:1         A[Del(C):R2]T  A[Del(C):R2]T
      19:     1      0        DEL:C:1:3     A[Del(C):R(4,5)]A  A[Del(C):R4]A
      20:     1      0  DEL:repeats:2:0 Del(2,):U(1,2):R(2,4)     Del2:U1:R2
      21:     3      0  DEL:repeats:3:1       Del(3,):U(3):R2     Del3:U3:R2
      22:     1      0       INS:T:1:5+     G[Ins(T):R(5,6)]G  G[Ins(T):R5]G
      23:    12      1      DEL:MH:5+:1            Del(6,):M1     Del(7,):M1
      24:     6      3      DEL:MH:5+:3            Del(6,):M3        Del6:M3
      25:     1      0       INS:C:1:5+          Ins(C):R(7,)  A[Ins(C):R7]A
      26:     1      0        DEL:C:1:1         G[Del(C):R2]A  G[Del(C):R2]A
      27:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      28:     1      0        DEL:C:1:3     A[Del(C):R(4,5)]A  A[Del(C):R4]A
      29:     1      0        DEL:T:1:2      C[Del(T):R(8,)]G  C[Del(T):R3]G
      30:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      31:     1      0        DEL:T:1:2      A[Del(T):R(8,)]A  A[Del(T):R3]A
      32:     1      0        DEL:T:1:0      G[Del(T):R(8,)]A  G[Del(T):R1]A
      33:     1      0        INS:T:1:4     C[Ins(T):R(0,4)]G  C[Ins(T):R4]G
      34:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      35:     1      0       INS:T:1:5+     G[Ins(T):R(5,6)]G  G[Ins(T):R6]G
      36:     1      0        INS:T:1:3     G[Ins(T):R(0,4)]C  G[Ins(T):R3]C
      37:     1      0        DEL:C:1:0         G[Del(C):R1]A  G[Del(C):R1]A
      38:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      39:     1      0       INS:T:1:5+      T[Ins(T):R(9,)]A T[Ins(T):R9+]A
      40:     1      0        DEL:C:1:0         T[Del(C):R1]T  T[Del(C):R1]T
      41:     1      0        INS:T:1:4     G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      42:     1      0        DEL:C:1:0         T[Del(C):R1]T  T[Del(C):R1]T
      43:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      44:     7      2      DEL:MH:5+:2            Del(6,):M2     Del(7,):M2
      45:     1      0        DEL:C:1:3     T[Del(C):R(4,5)]T  T[Del(C):R4]T
      46:     1      0       INS:T:1:5+      G[Ins(T):R(9,)]C G[Ins(T):R9+]C
      47:     1      0        DEL:T:1:0      C[Del(T):R(8,)]G  C[Del(T):R1]G
      48:     1      0        DEL:T:1:0      C[Del(T):R(8,)]C  C[Del(T):R1]C
      49:     1      0        INS:T:1:1     G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      50:     1      0        INS:T:1:4     G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      51:     1      0        DEL:C:1:1         G[Del(C):R2]T  G[Del(C):R2]T
      52:     1      0        DEL:C:1:2       Del(C):R(1,5)]G  A[Del(C):R3]G
      53:     1      0        DEL:T:1:3      A[Del(T):R(8,)]G  A[Del(T):R4]G
      54:     1      0        DEL:C:1:0         T[Del(C):R1]T  T[Del(C):R1]T
      55:     1      0       INS:T:1:5+     A[Ins(T):R(7,8)]A  A[Ins(T):R7]A
      56:     2      0  DEL:repeats:2:1 Del(2,):U(1,2):R(2,4)     Del2:U2:R2
      57:     1      0        DEL:T:1:1      A[Del(T):R(8,)]G  A[Del(T):R2]G
      58:     1      0        INS:T:1:3     A[Ins(T):R(0,4)]G  A[Ins(T):R3]G
      59:     3      1  INS:repeats:3:0           Ins(2,4):R0      Ins(5,):M
      60:     1      0  DEL:repeats:2:0 Del(2,):U(1,2):R(2,4)     Del2:U1:R2
              U koh_mh        COSMIC_83                Koh_89        Koh_476
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

