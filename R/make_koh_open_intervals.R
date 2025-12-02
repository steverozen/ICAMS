# Edit koh89 categories to change e.t. R(5,9) to R(5,) -- make them right-open
# intervals

make_koh_open_intervals = function(xx) {
  tmp = xx$Koh89.annotate.class
  tmp = gsub("Ins(2,):R(5,9)", "Ins(2,):R(5,)", x = tmp, fixed = TRUE) # Ins(C):R(7,)
  tmp = gsub("Ins(C):R(7,9)", "Ins(C):R(7,)", x = tmp, fixed = TRUE)
  tmp = gsub("[Del(T):R(8,9)]", "[Del(T):R(8,)]", x = tmp, fixed = TRUE)
  tmp = gsub("[Ins(T):R(8,9)]", "[Ins(T):R(8,)]", x = tmp, fixed = TRUE)
  tmp = gsub(
    "Del(3,):U(3,):R(3,)",
    "Del(3,):U(3,):R(3,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(2,):U(1,2):R(5,9)",
    "Del(2,):U(1,2):R(5,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(3,):U(3,):R(3,9)",
    "Del(3,):U(3,):R(3,)",
    x = tmp,
    fixed = TRUE
  )
  tmp = gsub(
    "Del(2,8):U(1,2):R(2,4)",
    "Del(2,):U(1,2):R(2,4)",
    x = tmp,
    fixed = TRUE
  )

  tmp
}
