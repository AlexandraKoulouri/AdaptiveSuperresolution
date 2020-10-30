function [SLE,Flow] = SigSimularityEMD(Str_e,Str_a,xe,xa);  
%REQUIRE SORTING!!!



% Earth Mover's Distance METRIC
    %ESTIMATE THE distance 
    Dc = zeros(length(Str_e), length(Str_a));
    for i=1:size(xe,1);
        for j=1:size(xa,1);
            Dc(i,j) = sqrt((xe(i,1)-xa(j,1))^2+(xe(i,2)-xa(j,2))^2);
           % Dc(j,i) = Dc(i,j) ;
        end
    end


Dc = Dc';
%extra variables for EMD estimation
extra_mass_penalty=0;%-1;
flowType= 3;


Str_e(Str_e<10e-4*max(Str_e)) = 0; %Remove very small entries
Str_a(Str_a<10e-4*max(Str_a)) = 0; %Remove very small entries
%
%normalize strengths and estimate powers
nP_e = Str_e/sum(Str_e);
nP_a = Str_a/sum(Str_a); 

[SLE, Flow]= emd_hat_mex_nes(nP_a,nP_e,Dc,extra_mass_penalty,flowType);
%SLE = SLE*1000;
