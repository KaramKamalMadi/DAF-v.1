
%function Alamouti
clc
clear

N_symbols = 1000;
SNR = 0:3:12;
nt = 2  
nr = 2 
Power = 1 
M =4  % constellation size


Constellation_temp = qammod(0:(M-1),M)               % size M
S_norm = sqrt(sum(abs(Constellation_temp).^2)/4)

all_possible_combination = 0:(M-1) 
Constellation = qammod(0:(M-1),M)/S_norm       % % normalised constellation


BER  = zeros(1,length(SNR))
no_errors = zeros(1,length(SNR))


%-- Power constraint
 p_norm = Power/S_norm;


 
% SNR Loop
for ind_snr = 1:length(SNR)
     
    % M-QAM modulation
    input_symbols = floor(rand(2,N_symbols)*M) % Random sequence of transmitted symbols
    
    Seq = qammod(input_symbols, M)/S_norm     % Random sequence of modulated transmitted symbols
           
    %Noise Variance
    sigma = sqrt (Power/(10^(SNR(ind_snr)/10)))  %Variance       
           
   % Monte carlo simulation
   for ind_symb = 1:N_symbols 
      
       
           H = (randn(nr,nt) + j*randn(nr,nt))/sqrt(2)  % rayleigh Channel
         
          S1_S2 = Seq(:,ind_symb)  %% Transmited modulated sequence 2 symbole QAM [S1 ; S2]
           

            S_ALM_1 = [S1_S2(1) ; S1_S2(2) ] 
            S_ALM_2 = [-S1_S2(2)';  S1_S2(1)']
            
            % Transmisson
            N_1 = sigma*(randn(nr,1)+j*randn(nr,1))*1/sqrt(2)   
            Y_t1 = p_norm*H*S_ALM_1  + N_1  % y1(t1) y2(t1)
            
            N_2 = sigma*(randn(nr,1)+j*randn(nr,1))*1/sqrt(2)    
            Y_t2 = p_norm*H*S_ALM_2  + N_2 %y1(t2) y2(t2) 
              
        
             % Reception : application of alamouti combiner y_1_hat, y_2_hat
        
             y_1_h = H(1,1)'*Y_t1(1) + H(1,2)*Y_t2(1)'+ H(2,2)*Y_t2(2)' + H(2,1)'*Y_t1(2)
             
             y_2_h = H(1,2)'*Y_t1(1) - H(1,1)*Y_t2(1)' + H(2,2)'*Y_t1(2) -  H(2,1)*Y_t2(2)'
            
                         
             
              % Maximum-Likelihood detection:
              
              Constellation_test = Constellation 
              
              
              for n=1: M
                  Distl(n) = abs(y_1_h -  Constellation_test(n))
                  Dist2(n) = abs(y_2_h -  Constellation_test(n))
              end
               
               
              [Minl idx1] = min(Distl) 
              [Min2 idx2] = min(Dist2)
             
              
             % Calculation of error numbers: 
             
             sequence_test = [all_possible_combination(idx1) ; all_possible_combination(idx2)]   % Estimated symbols [S_1_hat S_2_hat]
             no_errors(ind_snr)  =  no_errors(ind_snr) + biterr(input_symbols(:,ind_symb),sequence_test,log2(M))
%  
%            
            
   end %symboles

   BER(ind_snr)  = no_errors(ind_snr)/(N_symbols*2*log2(M))
   
   
end % SNR


 BER



figure(1);
clf;%Clear current figure window
 %retains the current plot and certain axes properties so that subsequent graphing commands add to the existing graph.
  semilogy(SNR,BER,'r');%creates a plot using a base 10 logarithmic scale for the y-axis and a linear scale for the x-axis
hold on
 grid on ;

 
      
 


 
 
 
 
 
 
 
 
 
 
 
 
