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
              U koh_mh        COSMIC_83    pev_COSMIC_83
          <int>  <int>           <char>           <char>
       1:     1      0        DEL:C:1:1        DEL:C:1:1
       2:     1      0        DEL:T:1:0        DEL:T:1:0
       3:     1      0        INS:T:1:4        INS:T:1:4
       4:     1      0        DEL:C:1:1        DEL:C:1:1
       5:     1      0        DEL:T:1:1        DEL:T:1:1
       6:     1      0        DEL:T:1:0        DEL:T:1:0
       7:     1      0        INS:T:1:0        INS:T:1:0
       8:     1      0        INS:C:1:0        INS:C:1:0
       9:     2      0 INS:repeats:2:5+ INS:repeats:2:5+
      10:     1      0        DEL:C:1:3        DEL:C:1:3
      11:     1      0        DEL:T:1:2        DEL:T:1:2
      12:     1      0        DEL:C:1:1        DEL:C:1:1
      13:     1      0        INS:T:1:4        INS:T:1:4
      14:     9      2      DEL:MH:5+:2      DEL:MH:5+:2
      15:     1      0        DEL:T:1:2        DEL:T:1:2
      16:     4      3       DEL:MH:4:3       DEL:MH:4:3
      17:     1      0        DEL:T:1:3        DEL:T:1:3
      18:     1      0        DEL:C:1:1        DEL:C:1:1
      19:     1      0        DEL:C:1:3        DEL:C:1:3
      20:     1      0  DEL:repeats:2:0  DEL:repeats:2:0
      21:     3      0  DEL:repeats:3:1  DEL:repeats:3:1
      22:     1      0       INS:T:1:5+       INS:T:1:5+
      23:    12      1      DEL:MH:5+:1      DEL:MH:5+:1
      24:     6      3      DEL:MH:5+:3      DEL:MH:5+:3
      25:     1      0       INS:C:1:5+       INS:C:1:5+
      26:     1      0        DEL:C:1:1        DEL:C:1:1
      27:     1      0        INS:T:1:1        INS:T:1:1
      28:     1      0        DEL:C:1:3        DEL:C:1:3
      29:     1      0        DEL:T:1:2        DEL:T:1:2
      30:     1      0        INS:T:1:1        INS:T:1:1
      31:     1      0        DEL:T:1:2        DEL:T:1:2
      32:     1      0        DEL:T:1:0        DEL:T:1:0
      33:     1      0        INS:T:1:4        INS:T:1:4
      34:     1      0        INS:T:1:1        INS:T:1:1
      35:     1      0       INS:T:1:5+       INS:T:1:5+
      36:     1      0        INS:T:1:3        INS:T:1:3
      37:     1      0        DEL:C:1:0        DEL:C:1:0
      38:     1      0        INS:T:1:1        INS:T:1:1
      39:     1      0       INS:T:1:5+       INS:T:1:5+
      40:     1      0        DEL:C:1:0        DEL:C:1:0
      41:     1      0        INS:T:1:4        INS:T:1:4
      42:     1      0        DEL:C:1:0        DEL:C:1:0
      43:     1      0        INS:T:1:1        INS:T:1:1
      44:     7      2      DEL:MH:5+:2      DEL:MH:5+:2
      45:     1      0        DEL:C:1:3        DEL:C:1:3
      46:     1      0       INS:T:1:5+       INS:T:1:5+
      47:     1      0        DEL:T:1:0        DEL:T:1:0
      48:     1      0        DEL:T:1:0        DEL:T:1:0
      49:     1      0        INS:T:1:1        INS:T:1:1
      50:     1      0        INS:T:1:4        INS:T:1:4
      51:     1      0        DEL:C:1:1        DEL:C:1:1
      52:     1      0        DEL:C:1:2        DEL:C:1:2
      53:     1      0        DEL:T:1:3        DEL:T:1:3
      54:     1      0        DEL:C:1:0        DEL:C:1:0
      55:     1      0       INS:T:1:5+       INS:T:1:5+
      56:     2      0  DEL:repeats:2:1  DEL:repeats:2:1
      57:     1      0        DEL:T:1:1        DEL:T:1:1
      58:     1      0        INS:T:1:3        INS:T:1:3
      59:     3      1  INS:repeats:3:0  INS:repeats:3:0
      60:     1      0  DEL:repeats:2:0  DEL:repeats:2:0
              U koh_mh        COSMIC_83    pev_COSMIC_83

