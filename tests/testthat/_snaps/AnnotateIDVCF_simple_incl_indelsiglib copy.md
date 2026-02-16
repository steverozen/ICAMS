# AnnotateIDVCF_sample

    Code
      avcf1
    Output
              U     R                 Koh_89        Koh_476        COSMIC_83
          <int> <int>                 <char>         <char>           <char>
       1:     1     2           [Del(C):R2]A  G[Del(C):R2]A        DEL:C:1:1
       2:     1     1      A[Del(T):R(1,4)]G  A[Del(T):R1]G        DEL:T:1:0
       3:     1     4      C[Ins(T):R(0,4)]A  C[Ins(T):R4]A        INS:T:1:4
       4:     1     2           [Del(C):R2]T  A[Del(C):R2]T        DEL:C:1:1
       5:     1     2      A[Del(T):R(1,4)]A  A[Del(T):R2]A        DEL:T:1:1
       6:     1     1      A[Del(T):R(1,4)]C  A[Del(T):R1]C        DEL:T:1:0
       7:     1     0      G[Ins(T):R(0,4)]C  G[Ins(T):R0]C        INS:T:1:0
       8:     1     0          Ins(C):R(0,3)  G[Ins(C):R0]A        INS:C:1:0
       9:     2    11          Ins(2,):R(5,) Ins2:U2:R(5,9) INS:repeats:2:5+
      10:     1     4       [Del(C):R(4,5)]T  A[Del(C):R4]T        DEL:C:1:3
      11:     1     3      A[Del(T):R(1,4)]A  A[Del(T):R3]A        DEL:T:1:2
      12:     1     2           [Del(C):R2]T  T[Del(C):R2]T        DEL:C:1:1
      13:     1     4      G[Ins(T):R(0,4)]C  G[Ins(T):R4]C        INS:T:1:4
      14:     9     1             Del(6,):M2     Del(7,):M2      DEL:MH:5+:2
      15:     1     3      A[Del(T):R(1,4)]G  A[Del(T):R3]G        DEL:T:1:2
      16:     4     1        Del(4,5):M(3,4)        Del4:M3       DEL:MH:4:3
      17:     1     4      C[Del(T):R(1,4)]C  C[Del(T):R4]C        DEL:T:1:3
      18:     1     2           [Del(C):R2]T  A[Del(C):R2]T        DEL:C:1:1
      19:     1     4       [Del(C):R(4,5)]A  A[Del(C):R4]A        DEL:C:1:3
      20:     1     2            Del(2,4):R1     Del2:U1:R1  DEL:repeats:2:0
      21:     3     2       Del(3,):U(3,):R2     Del3:U3:R2  DEL:repeats:3:1
      22:     1     5      G[Ins(T):R(5,7)]G  G[Ins(T):R5]G       INS:T:1:5+
      23:    12     1             Del(6,):M1     Del(7,):M1      DEL:MH:5+:1
      24:     6     1             Del(6,):M3        Del6:M3      DEL:MH:5+:3
      25:     1     7           Ins(C):R(7,)  A[Ins(C):R7]A       INS:C:1:5+
      26:     1     2           [Del(C):R2]A  G[Del(C):R2]A        DEL:C:1:1
      27:     1     1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      28:     1     4       [Del(C):R(4,5)]A  A[Del(C):R4]A        DEL:C:1:3
      29:     1     3      C[Del(T):R(1,4)]G  C[Del(T):R3]G        DEL:T:1:2
      30:     1     1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      31:     1     3      A[Del(T):R(1,4)]A  A[Del(T):R3]A        DEL:T:1:2
      32:     1     1      G[Del(T):R(1,4)]A  G[Del(T):R1]A        DEL:T:1:0
      33:     1     4      C[Ins(T):R(0,4)]G  C[Ins(T):R4]G        INS:T:1:4
      34:     1     1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      35:     1     6      G[Ins(T):R(5,7)]G  G[Ins(T):R6]G       INS:T:1:5+
      36:     1     3      G[Ins(T):R(0,4)]C  G[Ins(T):R3]C        INS:T:1:3
      37:     1     1           [Del(C):R1]A  G[Del(C):R1]A        DEL:C:1:0
      38:     1     1      G[Ins(T):R(0,4)]G  G[Ins(T):R1]G        INS:T:1:1
      39:     1    20       T[Ins(T):R(8,)]A T[Ins(T):R20]A       INS:T:1:5+
      40:     1     1           [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      41:     1     4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A        INS:T:1:4
      42:     1     1           [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      43:     1     1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      44:     7     1             Del(6,):M2     Del(7,):M2      DEL:MH:5+:2
      45:     1     4       [Del(C):R(4,5)]T  T[Del(C):R4]T        DEL:C:1:3
      46:     1     9       G[Ins(T):R(8,)]C  G[Ins(T):R9]C       INS:T:1:5+
      47:     1     1      C[Del(T):R(1,4)]G  C[Del(T):R1]G        DEL:T:1:0
      48:     1     1      C[Del(T):R(1,4)]C  C[Del(T):R1]C        DEL:T:1:0
      49:     1     1      G[Ins(T):R(0,4)]A  G[Ins(T):R1]A        INS:T:1:1
      50:     1     4      G[Ins(T):R(0,4)]A  G[Ins(T):R4]A        INS:T:1:4
      51:     1     2           [Del(C):R2]T  G[Del(C):R2]T        DEL:C:1:1
      52:     1     3       [Del(C):R(1,5)]G  A[Del(C):R3]G        DEL:C:1:2
      53:     1     4      A[Del(T):R(1,4)]G  A[Del(T):R4]G        DEL:T:1:3
      54:     1     1           [Del(C):R1]T  T[Del(C):R1]T        DEL:C:1:0
      55:     1     7      A[Ins(T):R(5,7)]A  A[Ins(T):R7]A       INS:T:1:5+
      56:     2     2 Del(2,8):U(1,2):R(2,4)     Del2:U2:R2  DEL:repeats:2:1
      57:     1     2      A[Del(T):R(1,4)]G  A[Del(T):R2]G        DEL:T:1:1
      58:     1     3      A[Ins(T):R(0,4)]G  A[Ins(T):R3]G        INS:T:1:3
      59:     3     0            Ins(2,4):R0     Ins(2,4):M  INS:repeats:3:0
      60:     1     2            Del(2,4):R1     Del2:U1:R1  DEL:repeats:2:0
              U     R                 Koh_89        Koh_476        COSMIC_83
          <int> <int>                 <char>         <char>           <char>

---

    Code
      t(avcf1)
    Output
                [,1]            [,2]                [,3]               
      U         " 1"            " 1"                " 1"               
      R         " 2"            " 1"                " 4"               
      Koh_89    "[Del(C):R2]A"  "A[Del(T):R(1,4)]G" "C[Ins(T):R(0,4)]A"
      Koh_476   "G[Del(C):R2]A" "A[Del(T):R1]G"     "C[Ins(T):R4]A"    
      COSMIC_83 "DEL:C:1:1"     "DEL:T:1:0"         "INS:T:1:4"        
                [,4]            [,5]                [,6]               
      U         " 1"            " 1"                " 1"               
      R         " 2"            " 2"                " 1"               
      Koh_89    "[Del(C):R2]T"  "A[Del(T):R(1,4)]A" "A[Del(T):R(1,4)]C"
      Koh_476   "A[Del(C):R2]T" "A[Del(T):R2]A"     "A[Del(T):R1]C"    
      COSMIC_83 "DEL:C:1:1"     "DEL:T:1:1"         "DEL:T:1:0"        
                [,7]                [,8]            [,9]              
      U         " 1"                " 1"            " 2"              
      R         " 0"                " 0"            "11"              
      Koh_89    "G[Ins(T):R(0,4)]C" "Ins(C):R(0,3)" "Ins(2,):R(5,)"   
      Koh_476   "G[Ins(T):R0]C"     "G[Ins(C):R0]A" "Ins2:U2:R(5,9)"  
      COSMIC_83 "INS:T:1:0"         "INS:C:1:0"     "INS:repeats:2:5+"
                [,10]              [,11]               [,12]          
      U         " 1"               " 1"                " 1"           
      R         " 4"               " 3"                " 2"           
      Koh_89    "[Del(C):R(4,5)]T" "A[Del(T):R(1,4)]A" "[Del(C):R2]T" 
      Koh_476   "A[Del(C):R4]T"    "A[Del(T):R3]A"     "T[Del(C):R2]T"
      COSMIC_83 "DEL:C:1:3"        "DEL:T:1:2"         "DEL:C:1:1"    
                [,13]               [,14]         [,15]              
      U         " 1"                " 9"          " 1"               
      R         " 4"                " 1"          " 3"               
      Koh_89    "G[Ins(T):R(0,4)]C" "Del(6,):M2"  "A[Del(T):R(1,4)]G"
      Koh_476   "G[Ins(T):R4]C"     "Del(7,):M2"  "A[Del(T):R3]G"    
      COSMIC_83 "INS:T:1:4"         "DEL:MH:5+:2" "DEL:T:1:2"        
                [,16]             [,17]               [,18]          
      U         " 4"              " 1"                " 1"           
      R         " 1"              " 4"                " 2"           
      Koh_89    "Del(4,5):M(3,4)" "C[Del(T):R(1,4)]C" "[Del(C):R2]T" 
      Koh_476   "Del4:M3"         "C[Del(T):R4]C"     "A[Del(C):R2]T"
      COSMIC_83 "DEL:MH:4:3"      "DEL:T:1:3"         "DEL:C:1:1"    
                [,19]              [,20]             [,21]             
      U         " 1"               " 1"              " 3"              
      R         " 4"               " 2"              " 2"              
      Koh_89    "[Del(C):R(4,5)]A" "Del(2,4):R1"     "Del(3,):U(3,):R2"
      Koh_476   "A[Del(C):R4]A"    "Del2:U1:R1"      "Del3:U3:R2"      
      COSMIC_83 "DEL:C:1:3"        "DEL:repeats:2:0" "DEL:repeats:3:1" 
                [,22]               [,23]         [,24]         [,25]          
      U         " 1"                "12"          " 6"          " 1"           
      R         " 5"                " 1"          " 1"          " 7"           
      Koh_89    "G[Ins(T):R(5,7)]G" "Del(6,):M1"  "Del(6,):M3"  "Ins(C):R(7,)" 
      Koh_476   "G[Ins(T):R5]G"     "Del(7,):M1"  "Del6:M3"     "A[Ins(C):R7]A"
      COSMIC_83 "INS:T:1:5+"        "DEL:MH:5+:1" "DEL:MH:5+:3" "INS:C:1:5+"   
                [,26]           [,27]               [,28]             
      U         " 1"            " 1"                " 1"              
      R         " 2"            " 1"                " 4"              
      Koh_89    "[Del(C):R2]A"  "G[Ins(T):R(0,4)]G" "[Del(C):R(4,5)]A"
      Koh_476   "G[Del(C):R2]A" "G[Ins(T):R1]G"     "A[Del(C):R4]A"   
      COSMIC_83 "DEL:C:1:1"     "INS:T:1:1"         "DEL:C:1:3"       
                [,29]               [,30]               [,31]              
      U         " 1"                " 1"                " 1"               
      R         " 3"                " 1"                " 3"               
      Koh_89    "C[Del(T):R(1,4)]G" "G[Ins(T):R(0,4)]G" "A[Del(T):R(1,4)]A"
      Koh_476   "C[Del(T):R3]G"     "G[Ins(T):R1]G"     "A[Del(T):R3]A"    
      COSMIC_83 "DEL:T:1:2"         "INS:T:1:1"         "DEL:T:1:2"        
                [,32]               [,33]               [,34]              
      U         " 1"                " 1"                " 1"               
      R         " 1"                " 4"                " 1"               
      Koh_89    "G[Del(T):R(1,4)]A" "C[Ins(T):R(0,4)]G" "G[Ins(T):R(0,4)]A"
      Koh_476   "G[Del(T):R1]A"     "C[Ins(T):R4]G"     "G[Ins(T):R1]A"    
      COSMIC_83 "DEL:T:1:0"         "INS:T:1:4"         "INS:T:1:1"        
                [,35]               [,36]               [,37]          
      U         " 1"                " 1"                " 1"           
      R         " 6"                " 3"                " 1"           
      Koh_89    "G[Ins(T):R(5,7)]G" "G[Ins(T):R(0,4)]C" "[Del(C):R1]A" 
      Koh_476   "G[Ins(T):R6]G"     "G[Ins(T):R3]C"     "G[Del(C):R1]A"
      COSMIC_83 "INS:T:1:5+"        "INS:T:1:3"         "DEL:C:1:0"    
                [,38]               [,39]              [,40]          
      U         " 1"                " 1"               " 1"           
      R         " 1"                "20"               " 1"           
      Koh_89    "G[Ins(T):R(0,4)]G" "T[Ins(T):R(8,)]A" "[Del(C):R1]T" 
      Koh_476   "G[Ins(T):R1]G"     "T[Ins(T):R20]A"   "T[Del(C):R1]T"
      COSMIC_83 "INS:T:1:1"         "INS:T:1:5+"       "DEL:C:1:0"    
                [,41]               [,42]           [,43]               [,44]        
      U         " 1"                " 1"            " 1"                " 7"         
      R         " 4"                " 1"            " 1"                " 1"         
      Koh_89    "G[Ins(T):R(0,4)]A" "[Del(C):R1]T"  "G[Ins(T):R(0,4)]A" "Del(6,):M2" 
      Koh_476   "G[Ins(T):R4]A"     "T[Del(C):R1]T" "G[Ins(T):R1]A"     "Del(7,):M2" 
      COSMIC_83 "INS:T:1:4"         "DEL:C:1:0"     "INS:T:1:1"         "DEL:MH:5+:2"
                [,45]              [,46]              [,47]              
      U         " 1"               " 1"               " 1"               
      R         " 4"               " 9"               " 1"               
      Koh_89    "[Del(C):R(4,5)]T" "G[Ins(T):R(8,)]C" "C[Del(T):R(1,4)]G"
      Koh_476   "T[Del(C):R4]T"    "G[Ins(T):R9]C"    "C[Del(T):R1]G"    
      COSMIC_83 "DEL:C:1:3"        "INS:T:1:5+"       "DEL:T:1:0"        
                [,48]               [,49]               [,50]              
      U         " 1"                " 1"                " 1"               
      R         " 1"                " 1"                " 4"               
      Koh_89    "C[Del(T):R(1,4)]C" "G[Ins(T):R(0,4)]A" "G[Ins(T):R(0,4)]A"
      Koh_476   "C[Del(T):R1]C"     "G[Ins(T):R1]A"     "G[Ins(T):R4]A"    
      COSMIC_83 "DEL:T:1:0"         "INS:T:1:1"         "INS:T:1:4"        
                [,51]           [,52]              [,53]              
      U         " 1"            " 1"               " 1"               
      R         " 2"            " 3"               " 4"               
      Koh_89    "[Del(C):R2]T"  "[Del(C):R(1,5)]G" "A[Del(T):R(1,4)]G"
      Koh_476   "G[Del(C):R2]T" "A[Del(C):R3]G"    "A[Del(T):R4]G"    
      COSMIC_83 "DEL:C:1:1"     "DEL:C:1:2"        "DEL:T:1:3"        
                [,54]           [,55]               [,56]                   
      U         " 1"            " 1"                " 2"                    
      R         " 1"            " 7"                " 2"                    
      Koh_89    "[Del(C):R1]T"  "A[Ins(T):R(5,7)]A" "Del(2,8):U(1,2):R(2,4)"
      Koh_476   "T[Del(C):R1]T" "A[Ins(T):R7]A"     "Del2:U2:R2"            
      COSMIC_83 "DEL:C:1:0"     "INS:T:1:5+"        "DEL:repeats:2:1"       
                [,57]               [,58]               [,59]            
      U         " 1"                " 1"                " 3"             
      R         " 2"                " 3"                " 0"             
      Koh_89    "A[Del(T):R(1,4)]G" "A[Ins(T):R(0,4)]G" "Ins(2,4):R0"    
      Koh_476   "A[Del(T):R2]G"     "A[Ins(T):R3]G"     "Ins(2,4):M"     
      COSMIC_83 "DEL:T:1:1"         "INS:T:1:3"         "INS:repeats:3:0"
                [,60]            
      U         " 1"             
      R         " 2"             
      Koh_89    "Del(2,4):R1"    
      Koh_476   "Del2:U1:R1"     
      COSMIC_83 "DEL:repeats:2:0"

