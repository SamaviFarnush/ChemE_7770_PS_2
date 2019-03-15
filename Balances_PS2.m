%Samavi Farnush Bint E Naser
%CHEME 7770
%PS#2
%19 Feb 2019

function dxdt=Balances_PS2(t,x,S,A)

%general parameters
cell_volume=0.58*10^-15;                                        %L/cell
av_number=6.023*10^23;
%%assuming equal RNAP allocation/gene
RNAP_available=0.24*8000*10^6/(av_number*cell_volume);      %micro-M, BIND 101441
%%assuming equal gene allocation per mRNA type
gene_concentration_Gp=200*10^6/(av_number*cell_volume);         %micro-M 
%%assuming equal ribosome allocation per protein type
ribosome_available=0.8*45100*10^6/(av_number*cell_volume);  %micro-M, BIND 101441
dry_weight=280*10^-15;                                          %gDW/cell BIND 103904


L_TX=1000;                       %nt
L_TL=333;                        %aa

e_x=42;                          %nt/s BIND 108487
e_L=14;                          %aa/s BIND 108487

K_EX_avg=e_x/L_TX;                 %1/s
K_EL_avg=e_L/L_TL;                 %1/s

K_IX=60/42;                                   %McClure min^-1                                                    
saturation_constant=(80-50)/60*(38-10)*K_IX;  %McClure  micro-M
K_IL=60/15;


%%mRNA1
%control term
if t<=60
    I=10;                             %mM
else
    I=0;
end
KX1=0.3;                              %mM
n=1.5;
W1X=0.0;                              %no residual P1
W1I=100;
fI=(I^n)/(KX1^n+I^n);             %non-dimensional
u1=(W1X+W1I*fI)/(1+W1X+W1I*fI);   %control term

%Specific rate of transcription
Lj_1=1200;                        %nt
KE_1=K_EX_avg*L_TX/Lj_1*60;       %1/min
tau_1=KE_1/K_IX;                               
r_x_1=KE_1*RNAP_available*gene_concentration_Gp/(saturation_constant*tau_1+(tau_1+1)*gene_concentration_Gp); %micro-M/min
T_X1=r_x_1*u1*cell_volume/dry_weight; %umol/gDW/min

%%protein1
%Specific rate of translation
Ll_1=1200/3;                        %aa
KE_3=K_EL_avg*L_TL/Ll_1*60;         %1/min
tau_3=KE_3/K_IL;                               
r_x_3=KE_3*ribosome_available*x(1)/(saturation_constant*cell_volume/dry_weight*tau_3+(tau_3+1)*x(1)); %u-M/min
T_L1=r_x_3*cell_volume/dry_weight; %micromol/gDW/min

%%mRNA2
%control term
KX2=1;    %umol/gDW
n=1.5;
W2X=0.0;                                 %broken circuit residue=0
W12=1;
W32=100;
f12=(x(4)^n)/(KX2^n+x(4)^n);             %non-dimensional
KX2=10; 
f32=(x(6)^n)/(KX2^n+x(6)^n);             %non-dimensional
u2=(W2X+W12*f12+W32*f32)/(1+W2X+W12*f12+W32*f32);   %control term

%Specific rate of transcription
Lj_2=2400;                        %nt
KE_2=K_EX_avg*L_TX/Lj_2*60;       %1/min
tau_2=KE_2/K_IX;                               
r_x_2=KE_2*RNAP_available*gene_concentration_Gp/(saturation_constant*tau_2+(tau_2+1)*gene_concentration_Gp); %micro-M/min
T_X2=r_x_2*u2*cell_volume/dry_weight; %umol/gDW/min

%%protein2
%Specific rate of translation
Ll_2=2400/3;                        %aa
KE_5=K_EL_avg*L_TL/Ll_2*60;         %1/min
tau_5=KE_5/K_IL;                               
r_x_5=KE_5*ribosome_available*x(2)/(saturation_constant*cell_volume/dry_weight*tau_5+(tau_5+1)*x(2)); %u-M/min
T_L2=r_x_5*cell_volume/dry_weight; %umol/gDW/min

%%mRNA3
%control term
KX3=1;                          %umol/gDW
n=1.5;
W3X=0.0;
W13=5;
W23=100;
f13=(x(4)^n)/(KX3^n+x(4)^n);                                  %non-dimensional
KX3=10; 
f23=(x(5)^n)/(KX2^n+x(5)^n);
u3=(W3X+W13*f13+W23*f23)/(1+W3X+W13*f13+W23*f23);             %control term

%Specific rate of transcription
Lj_3=600;                        %nt
KE_3=K_EX_avg*L_TX/Lj_3*60;      %1/min
tau_3=KE_3/K_IX;                               
r_x_3=KE_3*RNAP_available*gene_concentration_Gp/(saturation_constant*tau_3+(tau_3+1)*gene_concentration_Gp); %micro-M/min
T_X3=r_x_3*u3*cell_volume/dry_weight; %micromol/gDW/min

%%protein3
%Specific rate of translation
Ll_3=600/3;                        %aa
KE_6=K_EL_avg*L_TL/Ll_3*60;       %1/min
tau_6=KE_6/K_IL;                               
r_x_6=KE_6*ribosome_available*x(3)/(saturation_constant*cell_volume/dry_weight*tau_6+(tau_6+1)*x(3)); %micro-M/min
T_L3=r_x_6*cell_volume/dry_weight; %micromol/gDW/min


r=[T_X1; T_X2; T_X3; T_L1; T_L2; T_L3];

    
dxdt=A*x+S*r;

end