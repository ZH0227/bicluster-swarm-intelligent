function biClustResult = floc(k,data )
%FLexible Overlapped biClustering Algorithm
% Yang J, et al. 
% An improved biclustering method for analyzing gene expression profiles. 
% Int. J. Artif. Intell. T. 2005;14:771-790.
%
% Inputs
%   data        -   input data matrix
%   numBiclust  -   the number of biclusters searched
%   pRow        -   genes initial probability of membership to the 
%                   biclusters
%   pColumn     -   samples initial probability of membership to the 
%                   biclusters
%   resTh       -   residue threshold
%   minRow      -   minimal number of gene per bicluster
%   minCol      -   minimal number of conditions per bicluster
%   nIter       -   number of iterations
%   blocRow     -   a matrix indicating the directed initialisation for
%                   the genes
%   blocColumn  -   a matrix indicating the directed initialisation for
%                   the conditions
%
% Outputs:
%   biClustResult: A structure consisting of
%       RowxNum     - Logical Matrix which contains 1 in [i,j] if Row i is in Bicluster j
%       NumxCol     - Logical Matrix which contains 1 in [i,j] if Col j is in Bicluster i
%       ClusterNo   - Number of clusters
%       Clust       - Another structure array containing all clusters with
%                     their respective row and column indices.
%       MatRes      - A matrix describing the biclusters with the following
%                     columns:
%                     Residue, Volume, Genes(rows), Conditions(columns), 
%                     Row Variance
%
% See Also: RESIDU
%
% Author: Jayesh Kumar Gupta, 2013.
%
% Contact: Jayesh Kumar Gupta http://home.iitk.ac.in/~jayeshkg
%          Indian Institute of Technology, Kanpur, India
if data == 1
  disp('import data: leukemia_big.csv')
  file = importdata('../data/leukemia_big.csv');
  data = file.data; %size 7128x72
  minRow = 20;
  minCol = 6;
elseif data==2
  disp('import data: yeast.matrix.txt')
  data = importdata('../data/yeast.matrix.txt'); %size 2884x17
  minRow = 10;
  minCol = 4;
elseif data==3
  disp('import data: breast_Matlab.txt')
  data = importdata('../data/breast_Matlab.txt'); %size 13666x117
  minRow = 30;
  minCol = 10;
elseif data == 4
  disp('import data: ratStrain.txt')
  data = importdata('../data/ratStrain.txt'); %size  8799x122
  minRow = 30;
  minCol = 10;
elseif data == 5
  disp('import data: DLBCL.txt')
  data = importdata('../data/DLBCL.txt'); %size  12625x21
  minRow = 30;
  minCol = 10;
end

  numBiclust = k;

  pRow = 0.5;
  pColumn = pRow;

  resTh = [];


  nIter = 500;

  blocRow = [];

  blocColumn = [];



if isempty(resTh), resTh = residu(data)/10; end

vecData = data(:);
nbRows = size(data,1);
nbCols = size(data,2);

if ~isempty(blocRow) || ~isempty(blocColumn)
  numBiclust = max([size(blocRow,2) size(blocColumn,2) numBiclust]);
end

vecBicRow = zeros(nbRows, numBiclust);
vecBicCol = zeros(nbCols, numBiclust);

if isempty(blocRow)
  blocRow = zeros(nbRows, numBiclust);
else 
  vecBicRow(:,1:size(blocRow,2)) = blocRow; 
end
vecBlocRow = vecBicRow(:);
vecBicRow = vecBicRow(:);

if isempty(blocColumn)
  blocColumn = zeros(nbCols, numBiclust);
else 
  vecBicCol(:,1:size(blocColumn,2)) = blocColumn; 
end
vecBlocCol = vecBicCol(:);
vecBicCol = vecBicCol(:);

% Generate Random
rand1 = rand(1, numBiclust*nbRows);
rand2 = rand(1, numBiclust*nbCols);

vecBicRow(rand1 < pRow) = 1;
vecBicCol(rand2 < pColumn) = 1;

vecResvol_bic = zeros(numBiclust*4,1);

mfloc( vecData, nbRows, nbCols, vecBicRow, vecBicCol, ...
                  vecResvol_bic, resTh, numBiclust, minRow, minCol,...
                  nIter, vecBlocRow, vecBlocCol);
bicRow = reshape(vecBicRow, nbRows, numBiclust)'; %Row wise reshape
bicCol = reshape(vecBicCol, nbCols, numBiclust)'; %Row wise reshape

matResvol_bic = zeros(numBiclust, 5);
tempmat = reshape(vecResvol_bic, 4, numBiclust)'; %Row wise reshape
matResvol_bic(:,1:4) = tempmat(:,1:4);
for i=1:numBiclust
  v = var(data(find(bicRow(i,:)==1), find(bicCol(i,:)==1)),0,2);
  matResvol_bic(i,5) = mean(v(:));
end
% Columns of matResvol_bic
% --------------------------------------------------------
% | Residue | Volume | Genes | Conditions | Row Variance |
% --------------------------------------------------------

% Result
RowxNum = bicRow';
NumxCol = bicCol;
for i=1:numBiclust
 rows = find(RowxNum(:,i)>0);
 cols =find(NumxCol(i,:)>0);
 cluster(i) = struct('rows', rows, 'cols', cols);
end

% Returning result
biClustResult = struct('RowxNum',RowxNum, 'NumxCol',NumxCol, ...
                        'ClusterNo', numBiclust, 'Clust', cluster, 'MatRes', matResvol_bic);
end

