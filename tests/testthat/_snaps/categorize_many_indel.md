# categorize_many_indels

    Code
      yy
    Output
          ins_or_del    pre ins_or_del_seq unmutated_rep_count   post    mh
              <char> <char>         <char>               <num> <char> <int>
       1:          d      G              C                   2      A    NA
       2:          d      A              T                   1      G    NA
       3:          i      C              T                   4      A    NA
       4:          d      A              C                   2      T    NA
       5:          d      A              T                   2      A    NA
       6:          d      A              T                   1      C    NA
       7:          i      G              T                   0      C    NA
       8:          i      G              C                   0      A    NA
       9:          i      G             TA                  10      T    NA
      10:          d      A              C                   4      T    NA
      11:          d      A              T                   3      A    NA
      12:          d      T              C                   2      T    NA
      13:          i      G              T                   4      C    NA
      14:          d      A      TAGTTATAC                   1      T     2
      15:          d      A              T                   3      G    NA
      16:          d      G           CCCA                   1      C     3
      17:          d      C              T                   4      C    NA
      18:          d      A              C                   2      T    NA
      19:          d      A              C                   4      A    NA
      20:          d      T             AA                   1      G    NA
      21:          d      A            CTT                   2      T    NA
      22:          i      G              T                   5      G    NA
      23:          d      G   CCTGAGCTGTAC                   1      C     1
      24:          d      G         CTGTAC                   1      C     3
      25:          i      A              C                   7      A    NA
      26:          d      G              C                   2      A    NA
      27:          i      G              T                   1      G    NA
      28:          d      A              C                   4      A    NA
      29:          d      C              T                   3      G    NA
      30:          i      G              T                   1      G    NA
      31:          d      A              T                   3      A    NA
      32:          d      G              T                   1      A    NA
      33:          i      C              T                   4      G    NA
      34:          i      G              T                   1      A    NA
      35:          i      G              T                   6      G    NA
      36:          i      G              T                   3      C    NA
      37:          d      G              C                   1      A    NA
      38:          i      G              T                   1      G    NA
      39:          i      T              T                  20      A    NA
      40:          d      T              C                   1      T    NA
      41:          i      G              T                   4      A    NA
      42:          d      T              C                   1      T    NA
      43:          i      G              T                   1      A    NA
      44:          d      A        TGATTCT                   1      T     2
      45:          d      T              C                   4      T    NA
      46:          i      G              T                   9      C    NA
      47:          d      C              T                   1      G    NA
      48:          d      C              T                   1      C    NA
      49:          i      G              T                   1      A    NA
      50:          i      G              T                   4      A    NA
      51:          d      G              C                   2      T    NA
      52:          d      A              C                   3      G    NA
      53:          d      A              T                   4      G    NA
      54:          d      T              C                   1      T    NA
      55:          i      A              T                   7      A    NA
      56:          d      C             TG                   2      T    NA
      57:          d      A              T                   2      G    NA
      58:          i      A              T                   3      G    NA
      59:          i      C            TCT                   0      T     1
      60:          d      A             GG                   1      T    NA
          ins_or_del    pre ins_or_del_seq unmutated_rep_count   post    mh
                 COSMIC_83    pev_COSMIC_83
                    <char>           <char>
       1:        DEL:C:1:1        DEL:C:1:1
       2:        DEL:T:1:0        DEL:T:1:0
       3:        INS:T:1:4        INS:T:1:4
       4:        DEL:C:1:1        DEL:C:1:1
       5:        DEL:T:1:1        DEL:T:1:1
       6:        DEL:T:1:0        DEL:T:1:0
       7:        INS:T:1:0        INS:T:1:0
       8:        INS:C:1:0        INS:C:1:0
       9: INS:repeats:2:5+ INS:repeats:2:5+
      10:        DEL:C:1:3        DEL:C:1:3
      11:        DEL:T:1:2        DEL:T:1:2
      12:        DEL:C:1:1        DEL:C:1:1
      13:        INS:T:1:4        INS:T:1:4
      14:      DEL:MH:5+:2      DEL:MH:5+:2
      15:        DEL:T:1:2        DEL:T:1:2
      16:       DEL:MH:4:3       DEL:MH:4:3
      17:        DEL:T:1:3        DEL:T:1:3
      18:        DEL:C:1:1        DEL:C:1:1
      19:        DEL:C:1:3        DEL:C:1:3
      20:  DEL:repeats:2:0  DEL:repeats:2:0
      21:  DEL:repeats:3:1  DEL:repeats:3:1
      22:       INS:T:1:5+       INS:T:1:5+
      23:      DEL:MH:5+:1      DEL:MH:5+:1
      24:      DEL:MH:5+:3      DEL:MH:5+:3
      25:       INS:C:1:5+       INS:C:1:5+
      26:        DEL:C:1:1        DEL:C:1:1
      27:        INS:T:1:1        INS:T:1:1
      28:        DEL:C:1:3        DEL:C:1:3
      29:        DEL:T:1:2        DEL:T:1:2
      30:        INS:T:1:1        INS:T:1:1
      31:        DEL:T:1:2        DEL:T:1:2
      32:        DEL:T:1:0        DEL:T:1:0
      33:        INS:T:1:4        INS:T:1:4
      34:        INS:T:1:1        INS:T:1:1
      35:       INS:T:1:5+       INS:T:1:5+
      36:        INS:T:1:3        INS:T:1:3
      37:        DEL:C:1:0        DEL:C:1:0
      38:        INS:T:1:1        INS:T:1:1
      39:       INS:T:1:5+       INS:T:1:5+
      40:        DEL:C:1:0        DEL:C:1:0
      41:        INS:T:1:4        INS:T:1:4
      42:        DEL:C:1:0        DEL:C:1:0
      43:        INS:T:1:1        INS:T:1:1
      44:      DEL:MH:5+:2      DEL:MH:5+:2
      45:        DEL:C:1:3        DEL:C:1:3
      46:       INS:T:1:5+       INS:T:1:5+
      47:        DEL:T:1:0        DEL:T:1:0
      48:        DEL:T:1:0        DEL:T:1:0
      49:        INS:T:1:1        INS:T:1:1
      50:        INS:T:1:4        INS:T:1:4
      51:        DEL:C:1:1        DEL:C:1:1
      52:        DEL:C:1:2        DEL:C:1:2
      53:        DEL:T:1:3        DEL:T:1:3
      54:        DEL:C:1:0        DEL:C:1:0
      55:       INS:T:1:5+       INS:T:1:5+
      56:  DEL:repeats:2:1  DEL:repeats:2:1
      57:        DEL:T:1:1        DEL:T:1:1
      58:        INS:T:1:3        INS:T:1:3
      59:  INS:repeats:3:0  INS:repeats:3:0
      60:  DEL:repeats:2:0  DEL:repeats:2:0
                 COSMIC_83    pev_COSMIC_83

