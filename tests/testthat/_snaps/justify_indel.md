# justify_indel examples work correctly

    Code
      justify_indel("CAAAG", "CAAG", pos = 2, expected_delta = "A")
    Output
      $leftmost_pos
      [1] 2
      
      $del_str
      [1] "A"
      
      $del_str_plus
      [1] "CA"
      
      $edge_warning
      [1] FALSE
      
      $error
      NULL
      

---

    Code
      justify_indel("CAAAG", "CAAG", pos = 3, expected_delta = "A")
    Output
      $leftmost_pos
      [1] 2
      
      $del_str
      [1] "A"
      
      $del_str_plus
      [1] "CA"
      
      $edge_warning
      [1] FALSE
      
      $error
      NULL
      

---

    Code
      justify_indel("CACAG", "CAG", pos = 3, expected_delta = "CA")
    Output
      $leftmost_pos
      [1] 2
      
      $del_str
      [1] "AC"
      
      $del_str_plus
      [1] "CAC"
      
      $edge_warning
      [1] TRUE
      
      $error
      NULL
      

---

    Code
      justify_indel("TCACAG", "TCAG", pos = 4, expected_delta = "CA")
    Output
      $leftmost_pos
      [1] 2
      
      $del_str
      [1] "CA"
      
      $del_str_plus
      [1] "TCA"
      
      $edge_warning
      [1] FALSE
      
      $error
      NULL
      

---

    Code
      justify_indel("TCACAG", "TCAG", pos = 3, expected_delta = "AC")
    Output
      $leftmost_pos
      [1] 2
      
      $del_str
      [1] "CA"
      
      $del_str_plus
      [1] "TCA"
      
      $edge_warning
      [1] FALSE
      
      $error
      NULL
      

---

    {
      "error": "expected_delta = A != substr(long_str, pos, pos + del_len - 1)) = AC"
    }

