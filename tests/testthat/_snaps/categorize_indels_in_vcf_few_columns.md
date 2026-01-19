# categorize_indels_in_vcf

    Code
      uu
    Output
                 COSMIC_83                 Koh_89        Koh_476
                    <char>                 <char>         <char>
       1:        DEL:C:1:1           [Del(C):R2]A  G[Del(C):R2]A
       2:        DEL:T:1:0      A[Del(T):R(1,4)]G  A[Del(T):R1]G
       3:        INS:T:1:4      C[Ins(T):R(0,4)]A  C[Ins(T):R4]A
       4:        DEL:C:1:1           [Del(C):R2]T  A[Del(C):R2]T
       5:        DEL:T:1:1      A[Del(T):R(1,4)]A  A[Del(T):R2]A
       6:        DEL:T:1:0      A[Del(T):R(1,4)]C  A[Del(T):R1]C
       7:        INS:T:1:0      G[Ins(T):R(0,4)]C  G[Ins(T):R0]C
       8:        INS:C:1:0          Ins(C):R(0,3)  G[Ins(C):R0]A
       9: INS:repeats:2:5+          Ins(2,):R(5,) Ins2:U2:R(5,9)
      10:        DEL:C:1:3       [Del(C):R(4,5)]T  A[Del(C):R4]T
      11:        DEL:T:1:2      A[Del(T):R(1,4)]A  A[Del(T):R3]A
      12:        DEL:C:1:1           [Del(C):R2]T  T[Del(C):R2]T
      13:        INS:T:1:4      G[Ins(T):R(0,4)]C  G[Ins(T):R4]C
      14:      DEL:MH:5+:2             Del(6,):M2     Del(7,):M2
      15:        DEL:T:1:2      A[Del(T):R(1,4)]G  A[Del(T):R3]G
      16:       DEL:MH:4:3        Del(4,5):M(3,4)        Del4:M3
      17:        DEL:T:1:3      C[Del(T):R(1,4)]C  C[Del(T):R4]C
      18:        DEL:C:1:1           [Del(C):R2]T  A[Del(C):R2]T
      19:        DEL:C:1:3       [Del(C):R(4,5)]A  A[Del(C):R4]A
      20:  DEL:repeats:2:0            Del(2,4):R1     Del2:U1:R1
      21:  DEL:repeats:3:1       Del(3,):U(3,):R2     Del3:U3:R2
      22:       INS:T:1:5+      G[Ins(T):R(5,7)]G  G[Ins(T):R5]G
      23:      DEL:MH:5+:1             Del(6,):M1     Del(7,):M1
      24:      DEL:MH:5+:3             Del(6,):M3        Del6:M3
      25:       INS:C:1:5+           Ins(C):R(7,)  A[Ins(C):R7]A
      26:        DEL:C:1:1           [Del(C):R2]A  G[Del(C):R2]A
      27:        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      28:        DEL:C:1:3       [Del(C):R(4,5)]A  A[Del(C):R4]A
      29:        DEL:T:1:2      C[Del(T):R(1,4)]G  C[Del(T):R3]G
      30:        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      31:        DEL:T:1:2      A[Del(T):R(1,4)]A  A[Del(T):R3]A
      32:        DEL:T:1:0      G[Del(T):R(1,4)]A  G[Del(T):R1]A
      33:        INS:T:1:4      C[Ins(T):R(0,4)]G  C[Ins(T):R4]G
      34:        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      35:       INS:T:1:5+      G[Ins(T):R(5,7)]G  G[Ins(T):R6]G
      36:        INS:T:1:3      G[Ins(T):R(0,4)]C  G[Ins(T):R3]C
      37:        DEL:C:1:0           [Del(C):R1]A  G[Del(C):R1]A
      38:        INS:T:1:1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G
      39:       INS:T:1:5+       T[Ins(T):R(8,)]A T[Ins(T):R20]A
      40:        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      41:        INS:T:1:4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      42:        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      43:        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      44:      DEL:MH:5+:2             Del(6,):M2     Del(7,):M2
      45:        DEL:C:1:3       [Del(C):R(4,5)]T  T[Del(C):R4]T
      46:       INS:T:1:5+       G[Ins(T):R(8,)]C  G[Ins(T):R9]C
      47:        DEL:T:1:0      C[Del(T):R(1,4)]G  C[Del(T):R1]G
      48:        DEL:T:1:0      C[Del(T):R(1,4)]C  C[Del(T):R1]C
      49:        INS:T:1:1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A
      50:        INS:T:1:4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A
      51:        DEL:C:1:1           [Del(C):R2]T  G[Del(C):R2]T
      52:        DEL:C:1:2       [Del(C):R(1,5)]G  A[Del(C):R3]G
      53:        DEL:T:1:3      A[Del(T):R(1,4)]G  A[Del(T):R4]G
      54:        DEL:C:1:0           [Del(C):R1]T  T[Del(C):R1]T
      55:       INS:T:1:5+      A[Ins(T):R(5,7)]A  A[Ins(T):R7]A
      56:  DEL:repeats:2:1 Del(2,8):U(1,2):R(2,4)     Del2:U2:R2
      57:        DEL:T:1:1      A[Del(T):R(1,4)]G  A[Del(T):R2]G
      58:        INS:T:1:3      A[Ins(T):R(0,4)]G  A[Ins(T):R3]G
      59:  INS:repeats:3:0            Ins(2,4):R0     Ins(2,4):M
      60:  DEL:repeats:2:0            Del(2,4):R1     Del2:U1:R1
                 COSMIC_83                 Koh_89        Koh_476
                    <char>                 <char>         <char>

