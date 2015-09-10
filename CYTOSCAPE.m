

%% Initializing workspace

   close all force ;
   clear all ;
   clc ;
   addpath('FUNCTIONS') ;

%% Randomnes conrtol

   seed = 0 ;
   RandStream.setGlobalStream( RandStream('mt19937ar','Seed',seed) ) ;
   randn( 'seed' , seed ) ;
   rand( 'seed' , seed ) ;

%% Output log

   disp( ['==============================================='] ) ;
   disp( ['= Running GRN inference using the VBEM method ='] ) ;
   disp( ['==============================================='] ) ; 
   disp(char(10)) ; disp( ['( use [Ctrl]+[C] to abort the execution )'] ) ;
 
%% Loading the data

   load('data.mat') ; Y = Y.(':')(':') ;% Y(:,172:end) = [] ; samples(172:end) = [] ;
%   Y = Y - repmat(mean(Y,2),1,size(Y,2)) ;
%   Y = (Y-repmat(mean(Y,2),1,numel(samples)))
   tic ; disp(char(10)) ;
   disp( [' # Choosing dataset...'] ) ;
%    [ input_file,input_path ] = uigetfile( {'*.tsv','TSV (*.tsv)';'*.*','Any file (*.*)'},'MultiSelect','off' ) ;
%    file = strsplit( input_file,'.' ) ;
%    disp( file(1) ) ;
   G = size(Y,1) ;
   N = size(Y,2)-1 ;
   
%% Variables and parameters initialization

 % Performance
   model = 'AR1MA1' ;
   cl = 1e-10 ;
   threshold = 0.5 ;
 % Outputs
   X = nan(G) ;
   W = nan(G) ;
   GRN = cell(0,5) ;
 % Priors
   m_x = 0.5*ones(G,1) ;
   S_x = 0.25*eye(G) ;
   m_w = zeros(G,1) ;
   S_w = eye(G) ;
   a = 2 ;
   b = 1/a ;
 % Posteriors
   mu_x = cell(G,1) ;
   SIGMA_x = cell(G,1) ;
   mu_w = cell(G,1) ;
   SIGMA_w = cell(G,1) ;
   alpha = cell(G,1) ;
   beta = cell(G,1) ;
   
%% Hyperparameters learning

   toc ; disp(char(10)) ;
   for i = 1:G

      disp( strjoin([' - ',genes(i)],'') ) ;
      [ mu_x{i} , SIGMA_x{i} , mu_w{i} , SIGMA_w{i} , alpha{i} , beta{i} ] = HYPERPARAMETERS( model , cl , Y , i , m_x , S_x , m_w , S_w , a , b ) ;

   end%for
   
%% GRN inference

   toc ; disp(char(10)) ;
   for i = 1:G
      probability = POSTERIOR( mu_x{i},SIGMA_x{i} ) ;
      parents = find( probability >= threshold ) ;
      if ( sum(parents) > 0 )
         X(parents,i) = 1 ;
         W(parents,i) = mu_w{i}(parents) ;
         for j = 1:numel(parents)
            GRN = [ GRN ; { genes{parents(j)} '->' genes{i} mu_w{i}(parents(j)) probability(parents(j)) } ] ;
         end%for
      end%if
   end%for
   
%% Cytoscape output file

   output_file = 'GRN2.txt' ;
   fid = fopen(output_file,'w') ;
   fprintf( fid,'Parent -> Child Weight Probability\n') ;
   for k = 1:size(GRN,1)
      fprintf( fid,'%s %s %s %f %f\n',GRN{k,:} ) ;
   end%for
   fclose(fid) ;
  
