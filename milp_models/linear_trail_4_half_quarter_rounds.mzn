include "linear_model_chacha.mzn";
include "linear_trail_4_and_7.mzn";
int: word_size=32;
int: R=9;
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

constraint STATE[1, 0, 2, 31]=1;
constraint STATE[1, 1, 3, 31]=1;
constraint sum(STATE[1, 0..3, 0..3, 0..31]) = 2;
constraint STATE[0, 0, 2, 28]=1; % Added to do match with output of differential-linear distinguisher 
constraint sum(STATE[2, 0..3, 0..3, 0..31])=1;
%constraint STATE[8, 0, 0, 32-1-11]=1; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 0, 3, 32-1-15]=0; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 2, 3, 32-1-27]=0; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 2, 3, 32-1-30]=0; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 3, 3, 32-1-30]=1; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 2, 2, 32-1-25]=1; % Added to do match with output of differential-linear distinguisher 2
%constraint STATE[8, 1, 3, 32-1-7]=0; % Added to do match with output of differential-linear distinguisher 2

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



output [
    "round:"  ++ show(r) ++ "-->    " ++ show(concat([if show(STATE[r,i,j,k])="1" then 
      "$\\Delta^"++ show(r) ++ "_{" ++show(i*4+j)++ "," ++ show(31-k) ++ "},$ " 
    endif
    |i in 0..3, j in 0..3, k in 0..31]))++"\n" | r in 0..R
];




%output ["words: \n"];
%output [ if show(STATE[r,i,j,k])="1" then show(i*4+j) ++"\n" else "" endif |i in 0..3, j in 0..3, k in 0..31,r in 8..8];
%output ["bits: \n"];
%output [ if show(STATE[r,i,j,k])="1" then show(31-k) ++"\n" else "" endif |i in 0..3, j in 0..3, k in 0..31,r in 8..8];

%output [ if show(STATE[r,i,j,k])="1" then "$\\Delta^7"++ "_{" ++show(i*4+j)++ "," ++ show(31-k) ++ "},$ \n" else "" endif |i in 0..3, j in 0..3, k in 0..31,r in 8..8];
%output [ if show(STATE[r,i,j,k])="1" then "$\\Delta^7.5"++ "_{" ++show(i*4+j)++ "," ++ show(31-k) ++ "},$ \n" else "" endif |i in 0..3, j in 0..3, k in 0..31,r in 9..9];



output ["Probability weight at round 7.5: "++show(sum(d_list[0..R, 0..3, 0..1,0..n-1]))++"\n"];
%output ["Hamming Weight: "++show(sum(STATE[0, 0..3, 0..3,0..n-1]))];

%output [concat([concat(["constraint STATE[" ++ show(r)  ++ "," ++ show(i) ++ "," ++ show(j) ++ "," ++ show(k) ++ "]=" ++ show(STATE[r,i,j,k])++"; "|i in 0..3, j in 0..3, k in 0..31])|r in 0..R])];