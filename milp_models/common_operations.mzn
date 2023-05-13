
function array[0..word_size-1] of var -1.0..1.0: GetWordFloat(array [0..3, 0..3, 0..word_size-1] of var -1.0..1.0 : STATE_, int: Index) = 
let {
    int: i = Index div 4;
    int: j = Index mod 4;
} in
array1d(0..word_size-1, [STATE_[i,j,k] | k in 0..word_size-1]);

function array[int] of var 0..1: GetWordBin(array [int, int, int] of var 0..1 : STATE_, int: Index) = 
let {
    int: n = 32;
    int: i = Index div 4;
    int: j = Index mod 4;
} in
array1d(0..n-1, [STATE_[i,j,k] | k in 0..n-1]);

function array [int] of string : float_array_to_string_array(array [int] of var float: X) =
let {
    int: n = length(X);
} in 
array1d(0..n-1, [format(X[i])|i in 0..n-1]);

function array [int] of string : int_array_to_string_array(array [int] of var 0..1: X) =
let {
    int: n = length(X);
} in 
array1d(0..n-1, [format(X[i])|i in 0..n-1]);

function string : int_array_to_bin_str(array [int] of var 0..1: X) =
let {
    int: n = length(X);
} in 
join("", int_array_to_string_array(X));

% Left rotation of X by val positions
function array[int] of var 0..1: LRot(array[int] of var 0..1: X, int: val)=
let {
    int: n = length(X);
} in 
array1d(0..n-1, [X[(j+val) mod n] | j in 0..n-1]);

% Right rotation of X by val positions
function array[int] of var 0..1: RRot(array[int] of var 0..1: X, int: val)=
let {
    int: n = length(X);
} in 
array1d(0..n-1, [X[(n+j-val) mod n] | j in 0..n-1]);
