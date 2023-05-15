include "linear_model_chacha__linear_trail_3.mzn";
include "linear_trail_3_proof.mzn";
int: word_size=32;
int: R=7;
array [0..R, 0..3, 0..3, 0..(n-1)]  of var 0..1 : STATE;
array [0..R, 0..3, 0..12, 0..n-1]  of var 0..1 : Z; 
array [0..R, 0..3, 0..8, 0..n-1]  of var 0..1 : dummy;
array [0..R, 0..3, 0..1, 0..n] of var 0..1 : d_list; 
array [0..3, 0..3] of 0..15: COLUMNS=array2d(0..3, 0..3, 
[0, 4, 8, 12,
 1, 5, 9, 13,
 2, 6, 10, 14,
 3,7,11,15
 ]
 );
 
array [0..3, 0..3] of 0..63: DIAGONALS=array2d(0..3, 0..3,
[
 0, 5, 10, 15,
 1, 6, 11, 12,
 2, 7, 8, 13,
 3, 4, 9, 14
]);


constraint sum(STATE[0, 0..3, 0..3, 0..31]) = 1;

constraint LinearChaCha(
  STATE,
  Z,
  dummy,
  d_list,
  R,
  1
) = true;

solve:: int_search(d_list, smallest, indomain_min, complete)
minimize sum(d_list[0..R, 0..3, 0..1,0..n-1]);
output [ if show(STATE[r,i,j,k])="1" then "$\\Delta^0"++ "_{" ++show(i*4+j)++ "," ++ show(31-k) ++ "},$ \n" else "" endif |i in 0..3, j in 0..3, k in 0..31,r in 0..0];


output ["Correlation: "++show(sum(d_list[0..R, 0..3, 0..1,0..n-1]))++"\n"];


output [
    "round:"  ++ show(r) ++ "  -->    " ++ show(concat([if show(STATE[r,i,j,k])="1" then 
      "$\\Delta^"++ show(r) ++ "_{" ++show(i*4+j)++ "," ++ show(31-k) ++ "},$ " 
    endif
    |i in 0..3, j in 0..3, k in 0..31]))++"\n\n" | r in 0..R
];
