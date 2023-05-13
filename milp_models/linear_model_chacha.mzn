include "common_operations.mzn";
include "linear_operations.mzn";



int :n = 32;

predicate LinearChaCha(
  array [int, 0..3, 0..3, 0..(n-1)]  of var 0..1 : STATE_,
  array [int, 0..3, 0..12, 0..n-1]  of var 0..1 : Z_, 
  array [int, 0..3, 0..8, 0..n-1]  of var 0..1 : dummy_,
  array [int, 0..3, 0..1, 0..n] of var 0..1 : d_list_,
  int : R_,
  int : last_round
) =
  forall (i in 1..R_) (
   Round(
     STATE_[i-1, 0..3, 0..3, 0..n-1], 
     d_list_[i-1, 0..3, 0..1, 0..n], 
     STATE_[i, 0..3, 0..3, 0..n-1], 
     dummy_[i-1, 0..3, 0..8, 0..n-1], 
     Z_[i-1, 0..3, 0..12, 0..n-1],
     i-1,
     last_round
    )
   );
 
%constraint sum(STATE[0, 0..3, 0..3, 0..31])<32;
%constraint sum(STATE[3, 0..3, 0..3, 0..31])=1;


predicate Round(
  array [0..3, 0..3, 0..n-1] of var 0..1 : X, 
  array [0..3, 0..1, 0..n] of var 0..1 : p, 
  array [0..3, 0..3, 0..n-1] of var 0..1 : Y,  
  array [0..3, 0..8, 0..n-1] of var 0..1 : dummy_branch, 
  array [0..3, 0..12, 0..n-1] of var 0..1 : Z_,
  int: round,
  int: linearlastRound
) = 
let { 
var bool:HalfImpar1=forall(i in 0..3)(
qround_half1(
  GetWordBin(X, DIAGONALS[i,0]), GetWordBin(X, DIAGONALS[i,1]), GetWordBin(X, DIAGONALS[i,2]), GetWordBin(X, DIAGONALS[i,3]), 
  GetWordBin(Y, DIAGONALS[i,0]), GetWordBin(Y, DIAGONALS[i,1]), GetWordBin(Y, DIAGONALS[i,2]), GetWordBin(Y, DIAGONALS[i,3]),
  p[i,0..1,0..n], dummy_branch[i, 0..8, 0..n-1], Z_[i,0..12, 0..n-1] 
  )
);

var bool:HalfImpar2=forall(i in 0..3)(
qround_half2(
  GetWordBin(X, DIAGONALS[i,0]), GetWordBin(X, DIAGONALS[i,1]), GetWordBin(X, DIAGONALS[i,2]), GetWordBin(X, DIAGONALS[i,3]), 
  GetWordBin(Y, DIAGONALS[i,0]), GetWordBin(Y, DIAGONALS[i,1]), GetWordBin(Y, DIAGONALS[i,2]), GetWordBin(Y, DIAGONALS[i,3]),
  p[i,0..1,0..n], dummy_branch[i, 0..8, 0..n-1], Z_[i,0..12, 0..n-1] 
  )
);

var bool:HalfPar1=forall(i in 0..3)(
qround_half1(
  GetWordBin(X, COLUMNS[i,0]), GetWordBin(X, COLUMNS[i,1]), GetWordBin(X, COLUMNS[i,2]), GetWordBin(X, COLUMNS[i,3]), 
  GetWordBin(Y, COLUMNS[i,0]), GetWordBin(Y, COLUMNS[i,1]), GetWordBin(Y, COLUMNS[i,2]), GetWordBin(Y, COLUMNS[i,3]),
  p[i,0..1,0..n], dummy_branch[i, 0..8, 0..n-1], Z_[i,0..12, 0..n-1] 
  )
);

var bool:HalfPar2=forall(i in 0..3)(
qround_half2(
  GetWordBin(X, COLUMNS[i,0]), GetWordBin(X, COLUMNS[i,1]), GetWordBin(X, COLUMNS[i,2]), GetWordBin(X, COLUMNS[i,3]), 
  GetWordBin(Y, COLUMNS[i,0]), GetWordBin(Y, COLUMNS[i,1]), GetWordBin(Y, COLUMNS[i,2]), GetWordBin(Y, COLUMNS[i,3]),
  p[i,0..1,0..n], dummy_branch[i, 0..8, 0..n-1], Z_[i,0..12, 0..n-1] 
  )
);

} in 
 if (linearlastRound = 0) then 
   if(round = 0) then
     HalfImpar1
   elseif (((round-1) div 2) mod 2 = 0 /\ (round > 0)) then
      if((round-1) mod 2 = 0) then
       HalfPar2
      else
        HalfImpar1
      endif
   elseif(((round-1) div 2) mod 2 = 1 /\ (round > 0)) then
      if(((round-1) mod 2) = 0) then
        HalfImpar2
      else
        HalfImpar2
      endif
   else
     1=1
   endif
 else 
   if(round = 0) then
   HalfImpar1
   elseif (round=1) then
   HalfImpar2
      elseif (round=2) then
   HalfPar1
      elseif (round=3) then
   HalfPar2
      elseif (round=4) then
   HalfImpar1
      elseif (round=5) then
   HalfImpar2
      elseif (round=6) then
   HalfPar1
      elseif (round=7) then
   HalfPar2
      elseif (round=8) then
   HalfImpar1
   %if(round = 0) then
   %  HalfImpar1
   %elseif (((round-1) div 2) mod 2 = 0 /\ (round > 0)) then
   %   if(((round-1) mod 2) = 0) then
   %     HalfImpar2
   %   else
   %     HalfPar1
   %   endif
   %elseif(((round-1) div 2) mod 2 = 1 /\ (round > 0)) then
   %   if(((round-1) mod 2) = 0) then
   %     HalfPar2
   %   else
   %    HalfImpar1
   %   endif
   else
     1=1
   endif
endif;




predicate qround_half1(
        array [0..(n-1)] of var 0..1 : X0,
        array [0..(n-1)] of var 0..1 : X1,
        array [0..(n-1)] of var 0..1 : X2,
        array [0..(n-1)] of var 0..1 : X3,

        
        array [0..(n-1)] of var 0..1 : Y0,
        array [0..(n-1)] of var 0..1 : Y1,
        array [0..(n-1)] of var 0..1 : Y2,
        array [0..(n-1)] of var 0..1 : Y3,
        array [0..1, 0..(n)] of var 0..1 : p_,
        array [0..8, 0..n-1] of var 0..1 : dummy_branch_,
        array [0..12, 0..n-1] of var 0..1 : Zp,
        
       
) = 
linear_modular_addition_word(X0, Zp[6,0..n-1], Zp[1,0..n-1], p_[0,0..n]) /\
linear_modular_addition_word(X2, Zp[4,0..n-1], Zp[2,0..n-1], p_[1,0..n]) /\
three_way_fork_word(X1, Zp[6,0..n-1], RRot(Y1, 12) , dummy_branch_[1,0..n-1]) /\
three_way_fork_word(Zp[1,0..n-1], X3, Y0 , dummy_branch_[6,0..n-1]) /\
three_way_fork_word(LRot(X3, 16), Zp[4,0..n-1], Y3 , dummy_branch_[7,0..n-1]) /\
three_way_fork_word(Zp[2,0..n-1], RRot(Y1, 12), Y2 , dummy_branch_[8,0..n-1]);


predicate qround_half2(
        array [0..(n-1)] of var 0..1 : X0,
        array [0..(n-1)] of var 0..1 : X1,
        array [0..(n-1)] of var 0..1 : X2,
        array [0..(n-1)] of var 0..1 : X3,

        
        array [0..(n-1)] of var 0..1 : Y0,
        array [0..(n-1)] of var 0..1 : Y1,
        array [0..(n-1)] of var 0..1 : Y2,
        array [0..(n-1)] of var 0..1 : Y3,
        array [0..1, 0..(n)] of var 0..1 : p_,
        array [0..8, 0..n-1] of var 0..1 : dummy_branch_,
        array [0..12, 0..n-1] of var 0..1 : Zp,
        
       
) = 
linear_modular_addition_word(X0, Zp[3,0..n-1], Zp[9,0..n-1], p_[0,0..n]) /\
linear_modular_addition_word(X2, Zp[11,0..n-1], Zp[12,0..n-1], p_[1,0..n])/\
three_way_fork_word(X1, Zp[3,0..n-1], RRot(Y1, 7) , dummy_branch_[5,0..n-1]) /\
three_way_fork_word(Zp[9,0..n-1], X3, Y0 , dummy_branch_[6,0..n-1]) /\
three_way_fork_word(LRot(X3, 8), Zp[11,0..n-1], Y3 , dummy_branch_[7,0..n-1]) /\
three_way_fork_word(Zp[12,0..n-1], RRot(Y1, 7), Y2 , dummy_branch_[8,0..n-1]);


 
 
%predicate modular_addition_word(
% array [0..(n-1)] of var 0..1 : A, 
% array [0..(n-1)] of var 0..1 : B, 
% array [0..(n-1)] of var 0..1 : C, 
% array [0..n] of var 0..1 : di_list) = 
% forall (j in 0..n-1) (modular_addition_bit(di_list[j],C[j],A[j],B[j], di_list[j+1]))/\di_list[0] = 0;
%
%predicate three_way_fork_word(array [0..(n-1)] of var 0..1: X, array [0..(n-1)] of var 0..1: Y, array [0..(n-1)] of var 0..1: Z, array [0..(n-1)] of var 0..1: Dummy) = 
%forall (j in 0..n-1) (
%three_way_fork(X[j], Y[j], Z[j], Dummy[j])
%);
% 
%predicate three_way_fork(var 0..1:x, var 0..1:y, var 0..1:z, var 0..1:dummy) =
%  x + y + z - 2*dummy >= 0 /\
%  x + y + z <= 2 /\
%  dummy - x >= 0 /\
%  dummy - y >= 0 /\
%  dummy - z >= 0;
% 
% 
%predicate modular_addition_bit(var 0..1:di,var 0..1:a0,var 0..1:a1,var 0..1:a2, var 0..1:di1) =
% di - a0  - a1  +  a2  +  di1 >= 0 /\
% di + a0  + a1  -  a2  -  di1 >= 0 /\
% di + a0  - a1  -  a2  +  di1 >= 0 /\
% di - a0  + a1  -  a2  +  di1 >= 0 /\
% di + a0  - a1  +  a2  -  di1 >= 0 /\
% di - a0  + a1  +  a2  -  di1 >= 0 /\
% a0 - di  + a1  +  a2  +  di1 >= 0 /\
% di + a0  + a1  +  a2  +  di1 <= 4;
 
%function array[0..n-1] of var 0..1: LRot(array[0..n-1] of var 0..1: X, int: val)=
%array1d(0..n-1, [X[(j+val) mod n] | j in 0..n-1]);
%
%function array[0..n-1] of var 0..1: RRot(array[0..n-1] of var 0..1: X, int: val)=
%array1d(0..n-1, [X[(n+j-val) mod n] | j in 0..n-1]);

%solve:: int_search(d_list, smallest, indomain_min, complete)
%minimize sum(d_list[0..R, 0..3, 0..1,0..n-1]);

%output ["Round--------------->"++show(r)++"\n"++
%        "print('\\n'.join(textwrap.wrap("++
%        "'0x'+',0x'.join(textwrap.wrap('{0:#0{1}x}'.format(int('"++
%         concat([show(STATE[r,i,j,k])|i in 0..3, j in 0..3, k in 0..31])++
%         "',2), 130)[2:], 8))" ++
%         ",44)))\n"++ 
%        ", " ++ "\n"++ 
%        "State Weight"++":"++show(sum(STATE[r, 0..3, 0..3, 0..31]))++ "\n"|r in 0..R]; 
%output ["probability: "++show(sum(d_list[0..R,0..3, 0..1, 0..n-1]))++"\n"];
%%output ["probability round 1: "++show(sum(d_list[1,0..3, 0..1, 0..n-1]))++"\n"];
%%output ["probability round 1: "++show(sum(d_list[2,0..3, 0..1, 0..n-1]))++"\n"];
%%output ["probability round 1: "++show(sum(d_list[3,0..3, 0..1, 0..n-1]))++"\n"];
%
%output ["Probability: "++show(sum(d_list[0,0..3, 0..1, 0..n-1]))];
%output ["Probability: "++show(sum(d_list[1,0..3, 0..1, 0..n-1]))];
%output ["Probability: "++show(sum(d_list[2,0..3, 0..1, 0..n-1]))];
%output ["Probability: "++show(sum(d_list[3,0..3, 0..1, 0..n-1]))];
