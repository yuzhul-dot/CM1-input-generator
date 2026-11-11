% C.I.G (CM1 Input Generator)
% This code generates an input sounding file for CM1.
% ========================== Functions ================================== %
% rt.m ~ calculates water vapor mixing ratio.
% es.m ~ calculates equilibrium water vapor pressure.
% tskew.m ~ plots Skew-T plot. Modified by Matt Kumjian and Yuzhu Lin ba
% sed on the original code of MIT.
% ========================== References ================================= %
% [1] http://www2.mmm.ucar.edu/people/bryan/cm1/soundings/
% [2] Markowski, Paul, and Yvette Pamela Richardson. Mesoscale Meteorology
% in Midlatitudes. Wiley-Blackwell, 2010.
% [3] Warren, R.A., H. Richter, H.A. Ramsay, S.T. Siems, and M.J. Manton, 2
% 017: Impact of Variations in Upper-Level Shear on Simulated Supercells. M
% on.145, 2659?681, https://doi.org/10.1175/MWR-D-16-0412.1
% [4] Weisman, M. L., and R. Rotunno, 2000: The use of vertical wind shear
% versus helicity in interpreting supercell dynamics. J. Atmos. Sci., 57,
% 1452?472, doi:10.1175/1520-0469(2000)057,1452: TUOVWS.2.0.CO;2.
% ======================================================================= %
% Yuzhu Lin ~ yxl5930@psu.edu
% 11/04/2019 V1
% 05/20/2025 V2
%   /\ ___ /\
% \| = o _ o=| /3 meow meow meow meow
% ======================================================================= %

clear

%% variable allocation
% ======================================================================= %
%   These are the variables you would assign. If something breaks or the
%   skew-T looks off, check if the "levels" (Zml, LCL...) are set
%   appropriately.
%   If the codes indicate a "mismatch" it means your Tml is in conflict
%   with your LFC (e.g., Tml larger than parcel temp at Zml, but your LFC
%   is set above Zml. In this case, change either values, or Zml.
%   If the codes indicate EL cannot be achieved, this means your buoyancy
%   profile does not allow EL to be set that high in the stratosphere (
%   which is supposed to have constant T). You can either change your EL
%   height or your buoyancy profile parameters.
% ======================================================================= %

% ======================================================================= %
%   Variables without underline: variable value at a single level (e.g., Pg
%   represents atmospheric pressure at ground level, Tml represents the
%   temperature at mixed layer top). Variables with underline: columns of
%   variable values across levels, with increasing altitudes (e.g., T_ml
%   include temperature in the mixed layer, from bottom to top).
% ======================================================================= %

Zg=0;         %Ground Level height
Pg = 100000;  %Ground level Pressure
Tg = 25+273.15;     %Ground level temprature
%Hg (ground level Relative Humidity) is now calculated below so that LCL = Zml
Zml=1000;   %Mixed layer Height
lapse_ml = 7.5e-3; %Mixed layer lapse rate (K/m); g/Cpg (~9.67e-3) = well mixed (no CIN below Zml)
%Pml=88000;  %Mixed layer top Perssure -- now computed hydrostatically below (well-mixed theta_v)
Ztp=9500;	%Tropopause Height
Htp=0.90;   %Tropopause Relative Humity
Zst = 18000; %Domain top Height

% ======================================================================= %
%   Make VERY sure your LFC and EL are set appropriately in relations to
%   Zml and Ztp. For example, LFC can be higher or lower than Zml, but the
%   mixed layer lapse rate would be calculated differently.
% ======================================================================= %

%CIN = 150;   %CIN (calculated with virtual temperature) -- now an outcome of lapse_ml and LFC (case 2); still needed as input by cases 3 and 4
CAPE= 1500; %CAPE (calculated with virtual temperature)
Zb= 7000;   %Level of maximum buoyancy
CAPEblur = 5; % allowed maximum error for CAPE
LFC = 1500; %Level of free convection
EL = 9500; %Equilibrium level

ustrt=0; %u at 0m
vstrt=0; %v at 0m
% umid=10; %umax at 3km
% umax=30; %umax at 6km
vmax=7;  %umax at 6km
SRH = 100; %0-1 km SRH

%% constants

Rd=287.04;  %ideal gas constant for dry air
Cpd=1005;  %heat capacity for dry air
Cpv=1870;  %heat capacity for water vapor
Cpl=4218;  %heat capacity for liquid
lv=2.501e6;  %enthalpy for evaporazition
g=9.8;  %gravity constant
P0=100000;  %reference Pressure
epsi = 0.622; %Rd/Rv
lapsedry = 9.65e-3; %dry adiabatic lapse rate (reset to g/Cpg below for the surface parcel)

%% directly calculatable parameters
%   Surface RH that puts the LCL at Zml, so the parcel stays unsaturated (and
%   never buoyant) below Zml; CIN then comes from the mixed layer lapse rate
%   and the layer between Zml and LFC (see case 2). Saturation is placed half a grid step
%   (dzml/2 = 12.5 m) above Zml so the discrete LCL is the first level above Zml.
zsat = Zml + 12.5;
Hg = fzero(@(H) lcl_resid(H,Tg,Pg,zsat,Rd,Cpd,Cpv,g,P0,epsi), [0.05 1]);
esg=es(Tg); %ground equilibrium water vapor pressure
rtg =rt(Pg,esg,Hg); %ground water vapor mixing ratio
Cpg = Cpd*(1+rtg.*(Cpv./Cpd))./(1+rtg); %Cp at ground level
lapsedry = g./Cpg; %dry adiabat of the moist surface parcel: exactly neutral in the well-mixed (constant theta) layer
Rg = Rd.*(1+rtg./epsi)./(1+rtg); %R at ground level
Thetag = Tg*((P0./Pg)^(Rg./Cpg));   %PT at ground level
Tvg=Tg.*(1+rtg./epsi)./(1+rtg); %Ground Virtual Temperature
%theta at mixed layer top from the mixed layer lapse rate: dtheta/dz ~ (theta/T)(g/Cp - lapse_ml)
Thetaml = Thetag + (g/Cpg - lapse_ml)*Zml*Thetag/Tg;
rtml=rtg;  %moisture stays constant in the mixed layer
Thetavml = Thetag*(1+rtg/epsi)/(1+rtg);   %theta_v at the ground
betav = (Thetaml-Thetag)/Zml*(1+rtg/epsi)/(1+rtg);  %d(theta_v)/dz in the mixed layer
kap = Rd/Cpd;
%hydrostatic Exner function for linear theta_v: dpi/dz = -g/(Cpd*(Thetavml+betav*z))
if abs(betav) < 1e-12
    exn = @(z) (Pg/P0)^kap - g*z/(Cpd*Thetavml);
else
    exn = @(z) (Pg/P0)^kap - g/(Cpd*betav)*log(1+betav*z/Thetavml);
end
Pml = P0*exn(Zml)^(1/kap);  %hydrostatic mixed layer top pressure
Cpml = Cpd*(1+rtml.*(Cpv./Cpd))./(1+rtml);
Rml = Rd.*(1+rtml./epsi)./(1+rtml);
Tml=Thetaml/((P0./Pml)^(Rml./Cpml));
esml=es(Tml); %same stuff but at mixed layer top
Hml = rtml.*Pml/(esml*(0.622+ rtml));
Tvml=Tml.*(1+rtml./epsi)./(1+rtml);

%%%%%%%%%%%%%%%%%%%%%%%%%%%  vertical grids   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dzml = 25;                        %Mixed layer vertical Grids spacing
Z_ml=(0:dzml:Zml);                %Mixed Layer verticle Grids
mlgridsize=size(Z_ml);
%   You can absolutely change these to fit your needs. Just make sure they
%   start at the mixed layer top value you assigned earlier and ends at the
%   tropopause value you assigned earlier.
rt_ml=rtg;  %changed: well mixed;
P_ml = P0*exn(Z_ml).^(1/kap);  %hydrostatic, linear theta_v

dztp = 25;                        %Same stuff but for free troposphere
Z_tp=(Zml+dztp:dztp:Ztp);
tropogridsize=size(Z_tp);
%   Troposphere relative humidity profile- You can change this however you
%   want, just make sure H_tp starts at Hml and ends at Htp.
if Hml ==Htp
    H_tp =  zeros(1,tropogridsize(2))+Hml;
else
    H_tp = Hml+(Htp-Hml)/tropogridsize(2):(Htp-Hml)/tropogridsize(2):Htp;
end

dzst = 50;                       %Same stuff but for stratosphere
Z_st=(Ztp+dztp:dztp:Zst);
stgridsize = size(Z_st);
H_st=zeros(1,stgridsize(2))+Htp; %Assume H_st is constant

%% miscellaneous variables
%   Variables with _P denotes parcel properties
Z_P = transpose([Z_ml Z_tp Z_st]); %Parcel altitude
tempi = abs(Z_P-LFC);
iLFC = find(tempi == min(tempi),1); %LFC index
tempi = abs(Z_P-EL);
iEL = find(tempi == min(tempi),1);  %EL index
tempi = abs(Z_P-Zb);
iZb = find(tempi == min(tempi),1);  %Zb index

alli = mlgridsize(2) + tropogridsize(2) + stgridsize(2); %Domain size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parcel allocations %%%%%%%%%%%%%%%%%%%%%%%%%%
P_P = zeros(1,alli);                %Parcel pressure
rt_P = zeros(1,alli);               %Parcel mixing ratio
T_P = zeros(1,alli);                %Parcel Temperature
Tv_P = zeros(1,alli);               %Parcel Virtual Temperature
Theta_P = zeros(1,alli);            %Parcel Potential Temperature
es_P = zeros(1,alli);               %Parcel equilib. water vapor pressure
B_P = zeros(1,alli);                %Parcel buoyancy
R_P = zeros(1,alli);                %Parcel R

%%%%%%%%%%%%%%%%%%%%%%%%%%%% mixed layer allocations %%%%%%%%%%%%%%%%%%%%%%

Theta_ml=zeros(1,mlgridsize(2));
Tv_ml=zeros(1,mlgridsize(2));
T_ml=zeros(1,mlgridsize(2));
es_ml=zeros(1,mlgridsize(2));
H_ml=zeros(1,mlgridsize(2));
R_ml=zeros(1,mlgridsize(2));

P_tp=zeros(1,tropogridsize(2));
Theta_tp=zeros(1,tropogridsize(2));
Tv_tp=zeros(1,tropogridsize(2));
T_tp=zeros(1,tropogridsize(2));
es_tp=zeros(1,tropogridsize(2));
rt_tp=zeros(1,tropogridsize(2));
R_tp=zeros(1,tropogridsize(2));

P_st=zeros(1,stgridsize(2));
Theta_st=zeros(1,stgridsize(2));
Tv_st=zeros(1,stgridsize(2));
T_st=zeros(1,stgridsize(2));
es_st=zeros(1,stgridsize(2));
rt_st=zeros(1,stgridsize(2));
R_st=zeros(1,stgridsize(2));

%% LCL related calculations

%%%%%%%%%%%%%%%%%%%%%%%% Parcel peoperties below LCL %%%%%%%%%%%%%%%%%%%%%%

P_P(1: mlgridsize(2)) = P_ml; %P_P always equal to enviroment P

%   Initial parcel condition is the same as ground level atmospheric condition
T_P(1) = Tg;
Tv_P(1) = Tvg;
es_P(1) = esg;
Theta_P(1) = Thetag;

%   Check whether LCL is above or below Zml.
%   This is because pressure above Zml is unknown
rtes = 100;
i = 2;
while rtes > rtg
    if i> mlgridsize(2)
        P_P(i) = P_P(i-1)-50; %Placeholder for P above Zml, no real meaning
    end
    rt_P(i) = rtg; %rt_P equals to rtg below LCL
    CpP = Cpd*(1+rt_P(i).*(Cpv./Cpd))./(1+rt_P(i));
    R_P(i) = Rd.*(1+rt_P(i)./epsi)./(1+rt_P(i)); %Gas constant
    T_P(i) = T_P(i-1)-(Z_P(i)-Z_P(i-1)).*lapsedry; %Dry adiabatic lapse rate
    Tv_P(i) = T_P(i).*(1+rt_P(i)./epsi)./(1+rt_P(i)); %Dry adiabatic lapse rate
    es_P(i) = es(T_P(i)); %Equilib pressure
    rtes = rt(P_P(i),es_P(i),1); %Equilib mixed layer
    Theta_P(i) = T_P(i).*(P0./P_P(i)).^(R_P(i)./CpP); %Potential temperature
    i = i+1;
end

%%%%%%%%%%%%% Parcel properties below Zml for LCL <= Zml scenerio %%%%%%%%%

iLCL = i-1; %LCL index
LCL = Z_P(iLCL); %LCL height

%Parcel properties at LCL
if iLCL <= mlgridsize(2)
    PLCL_P = P_P(iLCL);
    TLCL_P = T_P(iLCL);
    TvLCL_P = Tv_P(iLCL);
    rtLCL_P = rt_P(iLCL);
    esLCL_P = es_P(iLCL);
    %Calculate the two variables below since they are the "anchor" for
    %parcel thermodynamic properties above LCL
    TsatLCL_P = 2840./(3.5.*log(TLCL_P)-log(esLCL_P/1000)-4.805)+55; %Saturation temperature
    ThetaepLCL_P = TLCL_P .*(P0./PLCL_P).^(0.2854.*(1-0.28.*rtLCL_P))...
        .*exp(rtLCL_P.*(1+0.81.*rtLCL_P).*(3376./TsatLCL_P-2.54)); %Pseudoequivalent PT
end

%Parcel properties between LCL and Zml
if iLCL<= mlgridsize(2)
    for i = iLCL:mlgridsize(2)
        %Solve for parcel T
        %Pseudoequibalent PT stays constant above LCL
        syms T
        Sml = ...
            vpasolve(ThetaepLCL_P == T .*...
            (P0./P_ml(i)).^(0.2854.*(1-0.28.*rt(P_ml(i),es(T),1)))...
            .*exp(rt(P_ml(i),es(T),1).*(1+0.81.*rt(P_ml(i),es(T),1)).*...
            (3376./(2840./(3.5.*log(T)-log(es(T)/1000)-4.805)+55)-2.54)),T);
        
        T_P(i) = double(Sml);
        es_P(i) = es(T_P(i));
        rt_P(i) = rt(P_P(i),es_P(i),1);
        
        Cp = Cpd*(1+rt_P(i).*(Cpv./Cpd))./(1+rt_P(i));
        Tv_P(i) = T_P(i).*(1+rt_P(i)./epsi)./(1+rt_P(i));
        R_P(i) = Rd.*(1+rt_P(i)./epsi)./(1+rt_P(i));
        Theta_P(i)=T_P(i)*(P0./P_P(i)).^(R_P(i)./Cp);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Below are buoyancy related calculations, divided into four cases      %
%   based on the relationship between LFC and Zml, as well as EL and      %
%   Ztp. This is set up for easier CAPE and CIN allocations.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% case 1
if LFC<=Zml && EL<=Ztp
    
    if Tv_P(mlgridsize(2))<Tvml
        'Tvml and LFC mismatch'
    end
    
    %%%%%%%%%%%%%%%%%%%% Enviroment properties below Zml %%%%%%%%%%%%%%%%%%
    
    % %solve lapse rate for mixed layer based on CIN, as LFC<Zml
    % %B(z) = g*(Tv_P - Tv)./Tv = a*z*(e^[b*(z-LFC)]-1)
    % %integration: a*(-z^2/2+(e^[b*(z-LFC)]*(b*z-1))/b^2)
    % %T_v = Tvg here.
    % syms a b
    % eqns = [a*Zml*(exp(b*(Zml-LFC))-1) == g.*(Tv_P(mlgridsize(2)) - Tvml)./Tvml,...
    %     -(a*LFC^2)/2 + (a*(-1+b*LFC + exp(-b*LFC)))*b^(-2)== -CIN];
    % S =vpasolve(eqns,[a b],[-5e-5 -1e-4]);
    % mla1 = double(S.a);
    % mlb1 = double(S.b);
    % 
    % if isempty(mla1) %if this formula does not work for our CIN, choose alternative
    %     ['alt buoyancy profile below Zml used']
    %     mla2 = CIN*6*LFC^(-3);
    %     B_P(1:iLFC) = mla2.*Z_P(1:iLFC).*(Z_P(1:iLFC)-LFC);
    %     B_P(iLFC:mlgridsize(2))= g.*(Tv_P(mlgridsize(2)) - Tvml)./Tvml./(Zml-LFC)*(Z_P(iLFC:mlgridsize(2))-LFC);
    % else
    %     B_P(1:mlgridsize(2)) = mla1.*Z_P(1:mlgridsize(2)).*(exp(mlb1.*(Z_P(1:mlgridsize(2))-LFC))-1);
    % end
    
    for i=1:mlgridsize(2)
        Theta_ml(i)=(Thetaml-Thetag).*Z_P(i)./Zml+Thetag;
        rt_ml(i) = rtml;
        R_ml(i) = Rd.*(1+rt_ml(i)./epsi)./(1+rt_ml(i));
        Cp = Cpd*(1+rt_ml(i).*(Cpv./Cpd))./(1+rt_ml(i));
        
        T_ml(i)=  Theta_ml(i)./((P0./P_ml(i)).^(R_ml(i)./Cp));
        Tv_ml(i)=  T_ml(i).*((1+rt_ml(i)./epsi)./(1+rt_ml(i)));
        es_ml(i)=es(T_ml(i));
        H_ml(i)=rt_ml(i)*P_ml(i)./(es(T_ml(i))*(epsi+rt_ml(i)));
        B_P(i) = (Tv_P(i) - Tv_ml(i))./Tv_ml(i).*g;
    end
    
    
    CINml = sum(B_P(1:iLFC).* dzml); %CIN below mixed layer top
    CAPEml = sum(B_P(iLFC:mlgridsize(2))).*dzml; %CAPE below mixed layer top
    
    PLFC_P = P_ml(iLFC);
    TLFC_P = T_P(iLFC);
    rtLFC_P = rt_P(iLFC);
    esLFC_P = es_P(iLFC);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Troposphere %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %solve buoyancy profile for troposphere based on CIN, as LFC<Zml
    %B(z) = tpa.*(H-z).*exp((z-H)./(H-Zb))-tpb.*(z-H).*(z-2.*Zb+H)
    %Integration of B(z) = CAPE-CAPEml
    Z_temp = Z_tp-Zml;
    Lscale = EL-Zml;
    Zbb = Zb-Zml;
    
    syms sca scb
    eqns = [B_P(mlgridsize(2)) == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale)))+2*scb*Zbb - scb*Lscale),...
        CAPE-CAPEml == (scb*Lscale^2*(-3*Zbb+2*Lscale))/3-...
        sca*(Zbb-Lscale)*(Zbb*(-1+exp(Lscale/(Zbb-Lscale)))+Lscale-2*exp(Lscale/(Zbb-Lscale))*Lscale)];
    Sab = vpasolve(eqns,[sca scb]);
    tpa = double(Sab.sca);
    tpb = double(Sab.scb);
    B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
        -tpb.*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
    
    %smooth out the B profile above EL- make it linear
    lapseEL = (B_P(iEL)- B_P(iEL-1))./dztp;
    B_P(iEL+1:tropogridsize(2)+mlgridsize(2)) = lapseEL.*(Z_P(iEL+1:tropogridsize(2)+mlgridsize(2))-EL);
    P_tp(1)= Pml.*exp( ...
        -g.*dztp./(Rd.*(1+rtml./epsi)./(1+rtml).* Tml));
    
    % Solve for troposphere properties, assume hydrostatic environment
    
    if LCL >Zml %LCL above Zml scenerio
        
        %%%% Below LCL %%%%
        rtes = 100; %equilib mixing ratio of parcel; temperary value here
        i = 1;
        while rtes > rtg
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            Z = Z_tp(i);
            
            %We are gonna solve the whole thing with a
            %7-unknown-variable-equation set :3
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),... %Hydrostatic: dp/dz = -rho*g
                rtP == rtg,... %Below LCL, parcel mixing ratio equalt to that at ground level (parcel origin)
                TP == Tml - lapsedry*(Z-Zml),... %Below LCL, parcel follows the adiabatic lapse rate
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),... %Parcel Tv calculared with
                rte == rt(P,es(Te),He),... %Enviromrntal mixing ratio
                Tve == TvP.*g./(g+Btemp),... %Environmental Tv calculated by parcel Tv and B
                Te ==  Tve./((1+rte./epsi)./(1+rte))]; %Environmental T calculated by environmental Tv
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            P_P(i+mlgridsize(2)) =  P_tp(i);
            
            % Other related calculations
            rtes = rt(P_tp(i),es(T_P(i+mlgridsize(2))),1);
            Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
            R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
            Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
            i = i+1;
        end
        
        %Parcel properties at LCL
        iLCL = i+mlgridsize(2)-1;
        LCL = Z_P(iLCL);
        PLCL_P = P_P(iLCL);
        TLCL_P = T_P(iLCL);
        TvLCL_P = Tv_P(iLCL);
        rtLCL_P = rt_P(iLCL);
        esLCL_P = es_P(iLCL);
        TsatLCL_P = 2840./(3.5.*log(TLCL_P)-log(esLCL_P/1000)-4.805)+55; %saturation temperature
        ThetaepLCL_P = TLCL_P .*(P0./PLCL_P).^(0.2854.*(1-0.28.*rtLCL_P))...
            .*exp(rtLCL_P.*(1+0.81.*rtLCL_P).*(3376./TsatLCL_P-2.54)); %pseudoequibalent PT
        
        %%%%% Above LCL %%%%%
        for i = iLCL-mlgridsize(2)+1:tropogridsize(2)
            
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            
            %We are gonna solve the whole thing with a
            %7-unknown-variable-equation set :3 only change is parcel rt
            %and parcel T calculation
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),... %Hydrostatic: dp/dz = -rho*g
                rtP == rt(P,es(TP),1),... %Parcel always at equilib mixing ratio above LCL
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54)),... %Parcel T calculated based on pseudoequibalent PT
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),... %Parcel Tv calculared with
                rte == rt(P,es(Te),He),... %Enviromrntal mixing ratio
                Tve == TvP.*g./(g+Btemp),... %Environmental Tv calculated by parcel Tv and B
                Te ==  Tve./((1+rte./epsi)./(1+rte))];%Environmental T calculated by environmental Tv
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            
            % Other related calculations
            P_P(i+mlgridsize(2)) =  P_tp(i);
            Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
            R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
            Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
        end
        
    else %LCL above Zml scenerio
        
        for i = 1:tropogridsize(2)
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            
            %We are gonna solve the whole thing with a
            %7-unknown-variable-equation set :3 only change is parcel rt
            %and parcel T calculation
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),... %Hydrostatic: dp/dz = -rho*g
                rtP == rt(P,es(TP),1),... %Parcel always at equilib mixing ratio above LCL
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54)),... %Parcel T calculated based on pseudoequibalent PT
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),... %Parcel Tv calculared with
                rte == rt(P,es(Te),He),... %Enviromrntal mixing ratio
                Tve == TvP.*g./(g+Btemp),... %Environmental Tv calculated by parcel Tv and B
                Te ==  Tve./((1+rte./epsi)./(1+rte))];%Environmental T calculated by environmental Tv
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            P_P(i+mlgridsize(2)) =  P_tp(i);
            
            % Other related calculations
            Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
            R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
            Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stratosphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    He = Htp;
    Te = T_tp(end); %isothermal
    P_st(1) = P_tp(end);
    
    for i = 1:stgridsize(2)
        
        Pund = P_st(max([i-1,1]));
        
        %same as troposphere, except that we don't need to calculate
        %enviromental temperature since stratosphere is isothermal
        syms P rte Tve TvP rtP TP
        eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
            rte == rt(P,es(Te),He),...
            Tve == Te.*((1+rte./epsi)./(1+rte)),...
            TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
            rtP == rt(P,es(TP),1),...
            ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
            .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
            -2.54))];
        
        Sall = vpasolve(eqns,[P rte Tve TvP rtP TP],[Pund;0.01;273.15;273.15;0.01;273.15]);
        
        T_st(i) = Te;
        P_st(i) = double(Sall.P);
        Tv_st(i) = double(Sall.Tve);
        rt_st(i) = double(Sall.rte);
        es_st(i)=es(T_st(i));
        
        rt_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.rtP);
        Tv_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TvP);
        T_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TP);
        P_P(i+tropogridsize(2)+ mlgridsize(2)) = P_st(i);
        
        % Other related calculations
        Cp = Cpd*(1+rt_st(i).*(Cpv./Cpd))./(1+rt_st(i));
        R_st(i) = Rd.*(1+rt_st(i)./epsi)./(1+rt_st(i));
        Theta_st(i)=T_st(i)*(P0./P_st(i)).^(R_st(i)./Cp);
    end
end

%% case 2

if LFC>Zml && EL<=Ztp
    if Tv_P(mlgridsize(2))>Tvml+0.01 %parcel is neutral at Zml in a well-mixed layer; allow round-off
        'Tvml and LFC mismatch'
    end
    
    CAPEml = 0; %CAPE below mixed layer top is 0 if LFC > Zml in our case
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Troposphere %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We calculate troposphere first, so that we know the CIN above Zml- so
    % that it can be deducted when calculating the B profile in the mixed
    % layer. Explanations for other calculations are commented in case 1.
    % solve buoyancy profile for troposphere based on CIN, as LFC<Zml
    % B(z) = tpa.*(H-z).*exp((z-H)./(H-Zb))-tpb.*(z-tpc).*(z-2*Zb+tpc)
    % Integration of B(z) = CAPE-CAPEml
    
    
    Z_temp = Z_tp-LFC;
    Lscale = EL-LFC;
    Zbb = Zb-LFC;
    syms sca scb
    eqns = [0 == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale))))-scb.*(-Lscale).*(-2*Zbb+Lscale),...
        CAPE-CAPEml == -scb*Lscale*(-3*Lscale^2 + 6*Lscale*Zbb + Lscale*(-3*Zbb+Lscale))/3-sca*(Zbb-Lscale)*(Zbb*(-1+exp(Lscale/(Zbb-Lscale)))+Lscale-2*exp(Lscale/(Zbb-Lscale))*Lscale)];
    Sab = vpasolve(eqns,[sca scb]);
    tpa = double(Sab.sca);
    tpb = double(Sab.scb);
    
    B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
        -tpb*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
    
    %smooth out the B profile above EL- make it linear
    lapseEL = (B_P(iEL)- B_P(iEL-1))./dztp;
    B_P(iEL+1:tropogridsize(2)+mlgridsize(2)) = lapseEL.*(Z_P(iEL+1:tropogridsize(2)+mlgridsize(2))-EL);
    P_tp(1)= Pml.*exp( ...
        -g.*dztp./(Rd.*(1+rtml./epsi)./(1+rtml).* Tml));
    lapseEL = (B_P(iEL)- B_P(iEL-1))./dztp;
    
    %Between Zml and LFC, B changes linearly from the mixed layer top
    %buoyancy B0 (continuous with the mixed layer) to 0 at LFC.
    %CIN is an outcome: CIN below Zml (set by lapse_ml) + (-B0)*(LFC-Zml)/2.
    B0 = g*(Tv_P(mlgridsize(2)) - Tvml)/Tvml;
    zcin = Z_P(mlgridsize(2):iLFC-1);
    B_P(mlgridsize(2):iLFC-1) = B0*(LFC-zcin)/(LFC-Zml);
    
    P_tp(1)= Pml.*exp( ...
        -g.*dztp./(Rd.*(1+rtml./epsi)./(1+rtml).* Tml));
    
    if LCL>Zml
        rtes = 100;
        i = 1;
        while rtes >= rtg
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            Z = Z_tp(i);
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                rte == rt(P,es(Te),He),...
                Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                Tve == TvP.*g./(g+Btemp),...
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                rtP == rtg,...
                TP == T_P(mlgridsize(2)) - lapsedry*(Z-Zml)];
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            P_P(i+mlgridsize(2)) =  P_tp(i);
            es_P(i+mlgridsize(2)) = es(T_P(i+mlgridsize(2)));
            rtes = rt(P_tp(i),es_P(i+mlgridsize(2)),1);
            i = i+1;
        end
        
        iLCL = i+mlgridsize(2)-1;
        LCL = Z_P(iLCL);
        PLCL_P = P_P(iLCL);
        TLCL_P = T_P(iLCL);
        TvLCL_P = Tv_P(iLCL);
        rtLCL_P = rt_P(iLCL);
        esLCL_P = es_P(iLCL);
        TsatLCL_P = 2840./(3.5.*log(TLCL_P)-log(esLCL_P/1000)-4.805)+55;
        ThetaepLCL_P = TLCL_P .*(P0./PLCL_P).^(0.2854.*(1-0.28.*rtLCL_P))...
            .*exp(rtLCL_P.*(1+0.81.*rtLCL_P).*(3376./TsatLCL_P-2.54));
        
        for i = iLCL-mlgridsize(2)+1:tropogridsize(2)
            
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                rte == rt(P,es(Te),He),...
                Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                Tve == TvP.*g./(g+Btemp),...
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                rtP == rt(P,es(TP),1),...
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54))];
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            P_P(i+mlgridsize(2)) =  P_tp(i);
            
            Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
            R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
            Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
        end
        
    else
        
        for i = 1:tropogridsize(2)
            
            Pund = P_tp(max([i-1,1]));
            Btemp = B_P(i+mlgridsize(2));
            He = H_tp(i);
            
            syms P rte Te Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                rte == rt(P,es(Te),He),...
                Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                Tve == TvP.*g./(g+Btemp),...
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                rtP == rt(P,es(TP),1),...
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54))];
            
            Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
            
            P_tp(i) = double(Sall.P);
            Tv_tp(i) = double(Sall.Tve);
            T_tp(i)= double(Sall.Te);
            rt_tp(i) = double(Sall.rte);
            es_tp(i)=es(T_tp(i));
            
            rt_P(i+mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
            T_P(i+mlgridsize(2)) = double(Sall.TP);
            P_P(i+mlgridsize(2)) =  P_tp(i);
            
            Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
            R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
            Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
        end
    end
    
    PLFC_P = P_tp(iLFC-mlgridsize(2));
    TLFC_P = T_P(iLFC);
    TvLFC_P = Tv_P(iLFC);
    rtLFC_P = rt_P(iLFC);
    esLFC_P = es_P(iLFC);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% above LCL, below Zml %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for i=1:mlgridsize(2)
        Theta_ml(i)=(Thetaml - Thetag).*Z_P(i)./Zml+Thetag;
        rt_ml(i) = rtml;
        R_ml(i) = Rd.*(1+rt_ml(i)./epsi)./(1+rt_ml(i));
        Cp = Cpd*(1+rt_ml(i).*(Cpv./Cpd))./(1+rt_ml(i));
        T_ml(i)=  Theta_ml(i)./((P0./P_ml(i)).^(R_ml(i)./Cp));
        Tv_ml(i)=  T_ml(i).*((1+rt_ml(i)./epsi)./(1+rt_ml(i)));
        es_ml(i)=es(T_ml(i));

        H_ml(i)=rt_ml(i)*P_ml(i)./(es(T_ml(i))*(epsi+rt_ml(i)));
        B_P(i) = (Tv_P(i) - Tv_ml(i))./Tv_ml(i).*g;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% above tropopause %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    He = Htp;
    Te = T_tp(end);
    P_st(1) = P_tp(end);
    for i = 1:stgridsize(2)
        
        Pund = P_st(max([i-1,1]));
        
        syms P rte Tve TvP rtP TP
        eqns = [P == Pund.*exp(-g.*dzst./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
            rte == rt(P,es(Te),He),...
            Tve == Te.*((1+rte./epsi)./(1+rte)),...
            TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
            rtP == rt(P,es(TP),1),...
            ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
            .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
            -2.54))]; %pseudoequibalent PT
        
        Sall = vpasolve(eqns,[P rte Tve TvP rtP TP],[Pund;0.01;273.15;273.15;0.01;273.15]);
        
        P_st(i) = double(Sall.P);
        Tv_st(i) = double(Sall.Tve);
        rt_st(i) = double(Sall.rte);
        T_st(i) = Te;
        es_st(i)=es(T_st(i));
        
        rt_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.rtP);
        Tv_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TvP);
        T_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TP);
        P_P(i+tropogridsize(2)+ mlgridsize(2)) = P_st(i);
        Cp = Cpd*(1+rt_st(i).*(Cpv./Cpd))./(1+rt_st(i));
        R_st(i) = Rd.*(1+rt_st(i)./epsi)./(1+rt_st(i));
        Theta_st(i)=T_st(i)*(P0./P_st(i)).^(R_st(i)./Cp);
    end
end

%% case 3

if LFC<=Zml && EL>Ztp
    
    if Tv_P(mlgridsize(2))<Tvml
        'Tvml and LFC mismatch'
    end
    
    %%%%%%%%%%%%%%%%%%%% Enviroment properties below Zml %%%%%%%%%%%%%%%%%%
    %see case 1; no change here
    syms a b
    eqns = [a*Zml*(exp(b*(Zml-LFC))-1) == g.*(Tv_P(mlgridsize(2)) - Tvml)./Tvml,...
        -(a*LFC^2)/2 + (a*(-1+b*LFC + exp(-b*LFC)))*b^(-2)== -CIN];
    S =vpasolve(eqns,[a b],[-5e-5 -1e-4]);
    mla1 = double(S.a);
    mlb1 = double(S.b);
    
    if isempty(mla1) %if this formula does not work for our CIN, choose alternative
        ['alt buoyancy profile below Zml used']
        mla2 = CIN*6*LFC^(-3);
        B_P(1:iLFC) = mla2.*Z_P(1:iLFC).*(Z_P(1:iLFC)-LFC);
        B_P(iLFC:mlgridsize(2))= g.*(Tv_P(mlgridsize(2)) - Tvml)./Tvml./(Zml-LFC)*(Z_P(iLFC:mlgridsize(2))-LFC);
    else
        B_P(1:mlgridsize(2)) = mla1.*Z_P(1:mlgridsize(2)).*(exp(mlb1.*(Z_P(1:mlgridsize(2))-LFC))-1);
    end
    
    for i=1:mlgridsize(2)
        Tv_ml(i) = Tv_P(i).*g./(B_P(i)+g);
        T_ml(i)=  Tv_ml(i)./((1+rt_ml(i)./epsi)./(1+rt_ml(i)));
        
        Cp = Cpd*(1+rt_ml(i).*(Cpv./Cpd))./(1+rt_ml(i));
        R_ml(i) = Rd.*(1+rt_ml(i)./epsi)./(1+rt_ml(i));
        Theta_ml(i)=T_ml(i)*(P0./P_ml(i)).^(R_ml(i)./Cp);
        es_ml(i)=es(T_ml(i));
        H_ml(i)=rt_ml(i)*P_ml(i)./(es(T_ml(i))*(epsi+rt_ml(i)));
        
    end
    
    CINml = sum(B_P(1:iLFC).* dzml);
    CAPEml = sum(B_P(iLFC:mlgridsize(2))).*dzml;
    
    PLFC_P = P_ml(iLFC);
    TLFC_P = T_P(iLFC);
    rtLFC_P = rt_P(iLFC);
    esLFC_P = es_P(iLFC);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% above Zml %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   Now, this part is a little tricky. Our stratosphere is isothermal,
    %   with T_st always equal to tropopause temperature. So there are
    %   three critirias we need to fullfill will out B profile:
    %       1. B = 0 at LFC, reaching maxima at Zb
    %       2. B reaches 0 at EL, which is in the isothermal stratosphere,
    %       whose environmental profile is entirely dependent on the
    %       condition at tropopause, which is calculated directly from
    %       troposphere buoyancy profile.
    %       3. CAPE in the troposphere and stratosphere add up together as our
    %       preassigned total CAPE- note that the stratosphere buoyancy
    %       profile is *decided* by the troposphere buoyancy profile (same
    %       reason as in critiria #2)
    %   These critirias demands initeration, to figure out how much
    %   buoyancy we put in troposphere. The method here is not exactly
    %   elegant- Hopefully I will change it in the future.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Troposphere %%%%%%%%%%%%%%%%%%%%%%%%%
    CAPEendloop = 0; %Triggered for the potential last CAPE iteration.
    CAPEendloopaccum = 0; %How many times potential last CAPE iteration has been triggered- for debugging.
    ELguesstrigger = 1; %Trigger for EL fine-tuning. When this is 1, we conduct iteration to reach EL.
    ELguesstriggeraccum = 0; %How many times EL tuning has been triggered- for debugging.
    Btpguess = 0; %Initual guess for buoyancy at Ztp
    CAPEstguess = Btpguess.*(EL-Ztp)./2; %Initual guess for stratosphere CAPE
    CAPEall = 0; %Calculated CAPE
    Z_temp = Z_tp-Zml;
    Ltp = Ztp-Zml;
    Zbb = Zb - Zml;
    %Lscale = EL+dzst-Zml; %Parameter changed in iterations that influences CAPE and EL
    ELBlapseguess = 1e-4; %Initial guess for buoyancy lapse rate in strat
    ELguess = Ztp+dzst; %Initial guess for calculated EL
    delEL = 100; %how much are we changing L per iteration step
    
    while abs(CAPEall-CAPE)>= CAPEblur || CAPEendloop ==1 %check if CAPE fullfill critiria 3
        %some EL iteration parameters
        Lscale = ELguess;
        ELCAPEstguess = CAPEstguess;
        
        ELCAPEstguessprev = ELCAPEstguess; %Document strat CAPE guess from the last iteration
        ELguessprev = ELguess; %Document EL guess from the last iteration
        Lprev = Lscale; %Document Lscale from the last iteration
        
        while ELguesstrigger<2
            ELCAPEstguessprev = ELCAPEstguess; %Document strat CAPE guess from the last iteration
            ELguessprev = ELguess; %Document EL guess from the last iteration
            Lprev = Lscale; %Document Lscale from the last iteration
            ELBlapseguessprev = ELBlapseguess;
            %calculate tropo buoyancy profile- same as case 1, but the
            %CAPE in strato is "guessed"
            
            Lscale = Lscale - (ELguess - EL);
            syms sca scb
            eqns = [B_P(mlgridsize(2)) == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale)))+2*scb*Zb - scb*Lscale),...
                CAPE-CAPEml - ELCAPEstguess == -(sca*exp((-Ltp+Lscale)/(Zbb-Lscale))*(Zbb-Lscale)*(Zbb*(-1+exp(Ltp/(Zbb-Lscale)))-Ltp-2*(-1+exp(Ltp/(Zbb-Lscale)))*Lscale))-(scb*Ltp*(-3*Zbb*Ltp+Ltp^2+6*Zbb*Lscale-3*Lscale^2))/3];
            Sab = vpasolve(eqns,[sca scb]);
            tpa = double(Sab.sca);
            tpb = double(Sab.scb);
            B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = ...
                tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
                -tpb.*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
            B_tp = B_P(tropogridsize(2)+mlgridsize(2));
            
            %Guess EL base on tropopause buoyancy and guessed buoyancy
            %lapse rate
            
%             ELBlapseguess = B_tp/(ELguess - Ztp);
%             ELguess = (ELCAPEstguess*2/ELBlapseguess)^0.5+Ztp;

            ELguess = ELCAPEstguess*2/B_tp+Z_tp;
            ELBlapseguess = B_tp/(ELguess - Ztp);
            
            %If guessed EL has approached allocated EL as much as it
            %can, stop the EL guess process (revert to best solution)
            if abs(ELguessprev-EL)-abs(ELguess-EL)<1e-3 && abs(ELguess-EL) <dzst/2
                Lscale = Lprev;
                %calculate troposphere buoyancy profile as usual, CAPE in the mixed
                %layer and stratosphere is deducted
                syms sca scb
                eqns = [B_P(mlgridsize(2)) == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale)))+2*scb*Zb - scb*Lscale),...
                    CAPE-CAPEml - ELCAPEstguessprev==   -(sca*exp((-Ltp+Lscale)/(Zbb-Lscale))*(Zbb-Lscale)*(Zbb*(-1+exp(Ltp/(Zbb-Lscale)))-Ltp-2*(-1+exp(Ltp/(Zbb-Lscale)))*Lscale))-(scb*Ltp*(-3*Zbb*Ltp+Ltp^2+6*Zbb*Lscale-3*Lscale^2))/3];
                Sab = vpasolve(eqns,[sca scb]);
                tpa = double(Sab.sca);
                tpb = double(Sab.scb);
                B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = ...
                    tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
                    -tpb.*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
                ELguess = ELCAPEstguess*2/B_tp+Z_tp;
                ELguesstrigger = 2; %set EL tuning trigger to 0 as we tried our best already
            end
        end
        
        %CAPEstguess = ELCAPEstguess; %Pass our new guess for strato CAPE
        if ELguesstrigger >0
            ['EL adjusting... EL now at ', num2str(ELguess)] %Annouce it to the class
        end
        % Tropo bupyancy profile calculation over
        
        %calculate tropo environment profile, same as case 1
        P_tp(1)= Pml.*exp( ...
            -g.*dztp./(Rd.*(1+rtml./epsi)./(1+rtml).* Tml));
        CAPEtp = sum(B_P(mlgridsize(2)+1:tropogridsize(2)+ mlgridsize(2)).*dztp);
        
        if LCL >Zml
            rtes = 100;
            i = 1;
            while rtes > rtg
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                Z = Z_tp(i);
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rtg,...
                    TP == Tml - lapsedry*(Z-Zml)]; %pseudoequibalent PT
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                rtes = rt(P_tp(i),es(T_P(i+mlgridsize(2))),1);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
                i = i+1;
            end
            
            iLCL = i+mlgridsize(2)-1;
            LCL = Z_P(iLCL);
            PLCL_P = P_P(iLCL);
            TLCL_P = T_P(iLCL);
            TvLCL_P = Tv_P(iLCL);
            rtLCL_P = rt_P(iLCL);
            esLCL_P = es_P(iLCL);
            TsatLCL_P = 2840./(3.5.*log(TLCL_P)-log(esLCL_P/1000)-4.805)+55;
            ThetaepLCL_P = TLCL_P .*(P0./PLCL_P).^(0.2854.*(1-0.28.*rtLCL_P))...
                .*exp(rtLCL_P.*(1+0.81.*rtLCL_P).*(3376./TsatLCL_P-2.54));
            
            for i = iLCL-mlgridsize(2)+1:tropogridsize(2)
                
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rt(P,es(TP),1),...
                    ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                    .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                    -2.54))]; %pseudoequibalent PT
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
            end
            
        else
            
            for i = 1:tropogridsize(2)
                
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rt(P,es(TP),1),...
                    ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                    .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                    -2.54))];
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% stratosphere %%%%%%%%%%%%%%%%%%%%%%%%
        %Same as case 1 and 2 here
        He = Htp;
        Te = T_tp(end);
        P_st(1) = P_tp(end);
        for i = 1:stgridsize(2)
            Pund = P_st(max([i-1,1]));
            syms P rte Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dzst./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                rte == rt(P,es(Te),He),...
                Tve == Te.*((1+rte./epsi)./(1+rte)),...
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                rtP == rt(P,es(TP),1),...
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54))]; %pseudoequibalent PT
            
            Sall = vpasolve(eqns,[P rte Tve TvP rtP TP],[Pund;0.01;273.15;273.15;0.01;273.15]);
            
            P_st(i) = double(Sall.P);
            Tv_st(i) = double(Sall.Tve);
            rt_st(i) = double(Sall.rte);
            T_st(i) = Te;
            es_st(i)=es(T_st(i));
            
            rt_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TvP);
            T_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TP);
            P_P(i+tropogridsize(2)+ mlgridsize(2)) = P_st(i);
            B_P(i+tropogridsize(2)+ mlgridsize(2)) = (Tv_P(i+tropogridsize(2)+ mlgridsize(2))-Tv_st(i))./Tv_st(i).*g;
            
            Cp = Cpd*(1+rt_st(i).*(Cpv./Cpd))./(1+rt_st(i));
            R_st(i) = Rd.*(1+rt_st(i)./epsi)./(1+rt_st(i));
            Theta_st(i)=T_st(i)*(P0./P_st(i)).^(R_st(i)./Cp);
            
            %Find EL
            if B_P(i+tropogridsize(2)+ mlgridsize(2)-1)>0 && B_P(i+tropogridsize(2)+ mlgridsize(2))<0
                iELnew = i+tropogridsize(2)+ mlgridsize(2);
                ELnew = Z_P(iELnew-1)*B_P(iELnew-1)/(B_P(iELnew-1)-B_P(iELnew))-Z_P(iELnew)*B_P(iELnew)/(B_P(iELnew-1)-B_P(iELnew));
            end
        end
        
        %Update some stuff
        
        CAPEst = sum(B_P(tropogridsize(2)+ mlgridsize(2)+1:iELnew).*dzst); %Actual calculated strato CAPE
        CAPEall = CAPEml + CAPEtp + CAPEst; %Actual calculated CAPE
        ELguess = ELnew;
        
        %New guess for stratosphere CAPE
        if abs(CAPEstguess - CAPEst)>=1
            CAPEstguess = (CAPEstguess+CAPEst)./2;
        else %if the CAPE iteration is stuck, force a pertubation
            CAPEstguess = CAPEst+(CAPEall-CAPE)*2;
        end
        
        %start trigger when calculated CAPE approach CAPE
        if abs(CAPEall-CAPE) < 1.5*CAPEblur
            ELguesstrigger = 1;
            ELguesstriggeraccum = ELguesstriggeraccum + ELguesstrigger;
        else
            ELguesstrigger = 0;
        end
        
        if abs(CAPEall-CAPE)<CAPEblur && CAPEendloopaccum<1
            CAPEendloop = 1;
            CAPEendloopaccum = CAPEendloopaccum + CAPEendloop;
        else
            CAPEendloop = 0;
            CAPEendloopaccum = 0;
        end
        
        ['still working on it... CAPE now at ',num2str(CAPEall)] %in case folks get anxious
        
    end
    
    
    %If EL simply can't be reached- usually because buoyancy simply cannot
    %be large enough at tropopause given the Zb and CAPE values, announce
    %it
    if abs(EL-ELnew) >= dzst
        ['EL unchievable - closest achievable EL chosen at ',num2str(ELnew),' m AGL']
    end
    
end
%% case 4

if LFC>Zml && EL>Ztp
    if Tv_P(mlgridsize(2))>=Tvml
        'Tvml and LFC mismatch'
    end
    
    CAPEml = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Troposphere %%%%%%%%%%%%%%%%%%%%%%%%%
    CAPEendloop = 0; %Triggered for the potential last CAPE iteration.
    CAPEendloopaccum = 0; %How many times potential last CAPE iteration has been triggered- for debugging.
    ELguesstrigger = 0; %Trigger for EL fine-tuning. When this is 1, we conduct iteration to reach EL.
    ELguesstriggeraccum = 0; %How many times EL tuning has been triggered- for debugging.
    Btpguess = 0; %Initual guess for buoyancy at Ztp
    CAPEstguess = Btpguess.*(EL-Ztp)./2; %Initual guess for stratosphere CAPE
    CAPEall = 0; %Calculated CAPE
    
    Z_temp = Z_tp-LFC;
    Ltp = Ztp - LFC;
    Zbb = Zb - LFC;
    ELBlapseguess = 1e-4; %Initial guess for buoyancy lapse rate in strat
    ELguess = Ztp+dzst; %Initial guess for calculated EL
    delEL = 100; %how much are we changing L per iteration step
    
    while abs(CAPEall-CAPE)>= CAPEblur  || CAPEendloop ==1 %check if CAPE fullfill critiria 3
        %some EL iteration parameters
        Lscale = ELguess;
        ELCAPEstguess = CAPEstguess;
        ELCAPEstguessprev = ELCAPEstguess; %Document strat CAPE guess from the last iteration
        ELguessprev = ELguess; %Document EL guess from the last iteration
        ELBlapseguessprev = ELBlapseguess;
        Lprev = Lscale; %Document Lscale from the last iteration
        
       Lscale = Lscale;
        
        while ELguesstrigger<2
            ELCAPEstguessprev = ELCAPEstguess; %Document strat CAPE guess from the last iteration
            ELguessprev = ELguess; %Document EL guess from the last iteration
            Lprev = Lscale; %Document Lscale from the last iteration
            ELBlapseguessprev = ELBlapseguess;
            Lscale = Lscale - (ELguess - EL);
            %calculate tropo buoyancy profile- same as case 1, but the
            %CAPE in strato is "guessed"
            syms sca scb
            eqns = [0 == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale))))-scb.*(-Lscale).*(-2*Zbb+Lscale),...
                CAPE-CAPEml - ELCAPEstguess==   -(sca*exp((-Ltp+Lscale)/(Zbb-Lscale))*(Zbb-Lscale)*(Zbb*(-1+exp(Ltp/(Zbb-Lscale)))-Ltp-2*(-1+exp(Ltp/(Zbb-Lscale)))*Lscale))-(scb*Ltp*(-3*Zbb*Ltp+Ltp^2+6*Zbb*Lscale-3*Lscale^2))/3];
            Sab = vpasolve(eqns,[sca scb]);
            tpa = double(Sab.sca);
            tpb = double(Sab.scb);
            
            B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = ...
                tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
                -tpb.*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
            B_tp = B_P(tropogridsize(2)+mlgridsize(2));
            %Guess EL base on tropopause buoyancy and guessed buoyancy
            %lapse rate
            ELBlapseguess = B_tp/(ELguess - Ztp);
            %ELguess = (ELCAPEstguess*2/ELBlapseguess)^0.5+Ztp;
            ELguess = ELCAPEstguess*2/B_tp+Ztp;
            %If guessed EL has approached allocated EL as much as it
            %can, stop the EL guess process (revert to best solution)
            if abs(ELguessprev-EL)-abs(ELguess-EL)<1e-3
                Lscale = Lprev;
                %calculate troposphere buoyancy profile as usual, CAPE in the mixed
                %layer and stratosphere is deducted
                syms sca scb
                eqns = [0 == -Lscale*(sca*(-exp(Lscale/(Zbb-Lscale))))-scb.*(-Lscale).*(-2*Zbb+Lscale),...
                    CAPE-CAPEml - ELCAPEstguessprev==   -(sca*exp((-Ltp+Lscale)/(Zbb-Lscale))*(Zbb-Lscale)*(Zbb*(-1+exp(Ltp/(Zbb-Lscale)))-Ltp-2*(-1+exp(Ltp/(Zbb-Lscale)))*Lscale))-(scb*Ltp*(-3*Zbb*Ltp+Ltp^2+6*Zbb*Lscale-3*Lscale^2))/3];
                Sab = vpasolve(eqns,[sca scb]);
                tpa = double(Sab.sca);
                tpb = double(Sab.scb);
                Z_temp = Z_tp-LFC;
                B_P(1+mlgridsize(2):tropogridsize(2)+mlgridsize(2)) = ...
                    tpa.*(Lscale-Z_temp).*exp((Z_temp-Lscale)./(Lscale-Zbb))...
                    -tpb.*(Z_temp-Lscale).*(Z_temp-2.*Zbb+Lscale);
                %ELguess = (ELCAPEstguessprev*2/ELBlapseguessprev)^0.5+Ztp;
                ELguess = ELCAPEstguess*2/B_tp+Ztp;
                ELguesstrigger = 2; %set EL tuning trigger to 0 as we tried our best already
            end
        end
        if ELguesstrigger >0
            ['EL adjusting... EL now at ', num2str(ELguess)] %Annouce it to the class
        end
        lapseLFC = -(g.*(Tv_P(mlgridsize(2)) - Tvml)/Tvml)./(LFC-Zml);
        B_P(mlgridsize(2):iLFC-1) = (g.*(Tv_P(mlgridsize(2)) - Tvml)/Tvml)+lapseLFC.*(Z_P(mlgridsize(2):iLFC-1)-Zml);
        %calculate tropo environment profile, same as case 1
        P_tp(1)= Pml.*exp( ...
            -g.*dztp./(Rd.*(1+rtml./epsi)./(1+rtml).* Tml));
        CAPEtp = sum(B_P(iLFC:tropogridsize(2)+ mlgridsize(2)).*dztp);
        
        if LCL >Zml
            rtes = 100;
            i = 1;
            while rtes > rtg
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                Z = Z_tp(i);
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rtg,...
                    TP == Tml - lapsedry*(Z-Zml)]; %pseudoequibalent PT
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                rtes = rt(P_tp(i),es(T_P(i+mlgridsize(2))),1);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
                i = i+1;
            end
            
            iLCL = i+mlgridsize(2)-1;
            LCL = Z_P(iLCL);
            PLCL_P = P_P(iLCL);
            TLCL_P = T_P(iLCL);
            TvLCL_P = Tv_P(iLCL);
            rtLCL_P = rt_P(iLCL);
            esLCL_P = es_P(iLCL);
            TsatLCL_P = 2840./(3.5.*log(TLCL_P)-log(esLCL_P/1000)-4.805)+55;
            ThetaepLCL_P = TLCL_P .*(P0./PLCL_P).^(0.2854.*(1-0.28.*rtLCL_P))...
                .*exp(rtLCL_P.*(1+0.81.*rtLCL_P).*(3376./TsatLCL_P-2.54));
            
            for i = iLCL-mlgridsize(2)+1:tropogridsize(2)
                
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rt(P,es(TP),1),...
                    ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                    .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                    -2.54))]; %pseudoequibalent PT
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
            end
            
        else
            
            for i = 1:tropogridsize(2)
                
                Pund = P_tp(max([i-1,1]));
                Btemp = B_P(i+mlgridsize(2));
                He = H_tp(i);
                
                syms P rte Te Tve TvP rtP TP
                eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                    rte == rt(P,es(Te),He),...
                    Te ==  Tve./((1+rte./epsi)./(1+rte)),...
                    Tve == TvP.*g./(g+Btemp),...
                    TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                    rtP == rt(P,es(TP),1),...
                    ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                    .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                    -2.54))];
                
                Sall = vpasolve(eqns,[P rte Te Tve TvP rtP TP],[Pund;0.01;273.15;273.15;273.15;0.01;273.15]);
                
                P_tp(i) = double(Sall.P);
                Tv_tp(i) = double(Sall.Tve);
                T_tp(i)= double(Sall.Te);
                rt_tp(i) = double(Sall.rte);
                es_tp(i)=es(T_tp(i));
                
                rt_P(i+mlgridsize(2)) = double(Sall.rtP);
                Tv_P(i+mlgridsize(2)) = double(Sall.TvP);
                T_P(i+mlgridsize(2)) = double(Sall.TP);
                P_P(i+mlgridsize(2)) =  P_tp(i);
                
                Cp = Cpd*(1+rt_tp(i).*(Cpv./Cpd))./(1+rt_tp(i));
                R_tp(i) = Rd.*(1+rt_tp(i)./epsi)./(1+rt_tp(i));
                Theta_tp(i)=T_tp(i)*(P0./P_tp(i)).^(R_tp(i)./Cp);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%% stratosphere %%%%%%%%%%%%%%%%%%%%%%%%
        %Same as case 1 and 2 here
        He = Htp;
        Te = T_tp(end);
        P_st(1) = P_tp(end);
        for i = 1:stgridsize(2)
            Pund = P_st(max([i-1,1]));
            syms P rte Tve TvP rtP TP
            eqns = [P == Pund.*exp(-g.*dztp./(Rd.*(1+rte./epsi)./(1+rte).* Te)),...
                rte == rt(P,es(Te),He),...
                Tve == Te.*((1+rte./epsi)./(1+rte)),...
                TvP == TP.*((1+rtP./epsi)./(1+rtP)),...
                rtP == rt(P,es(TP),1),...
                ThetaepLCL_P == TP .*(P0./P).^(0.2854.*(1-0.28.*rtP))...
                .*exp(rtP.*(1+0.81.*rtP).*(3376./(2840./(3.5.*log(TP)-log(es(TP)/1000)-4.805)+55)...
                -2.54))]; %pseudoequibalent PT
            
            Sall = vpasolve(eqns,[P rte Tve TvP rtP TP],[Pund;0.01;273.15;273.15;0.01;273.15]);
            
            P_st(i) = double(Sall.P);
            Tv_st(i) = double(Sall.Tve);
            rt_st(i) = double(Sall.rte);
            T_st(i) = Te;
            es_st(i)=es(T_st(i));
            
            rt_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.rtP);
            Tv_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TvP);
            T_P(i+tropogridsize(2)+ mlgridsize(2)) = double(Sall.TP);
            P_P(i+tropogridsize(2)+ mlgridsize(2)) = P_st(i);
            B_P(i+tropogridsize(2)+ mlgridsize(2)) = (Tv_P(i+tropogridsize(2)+ mlgridsize(2))-Tv_st(i))./Tv_st(i).*g;
            
            Cp = Cpd*(1+rt_st(i).*(Cpv./Cpd))./(1+rt_st(i));
            R_st(i) = Rd.*(1+rt_st(i)./epsi)./(1+rt_st(i));
            Theta_st(i)=T_st(i)*(P0./P_st(i)).^(R_st(i)./Cp);
            
            %Find EL
            if B_P(i+tropogridsize(2)+ mlgridsize(2)-1)>0 && B_P(i+tropogridsize(2)+ mlgridsize(2))<0
                iELnew = i+tropogridsize(2)+ mlgridsize(2);
                ELnew = Z_P(iELnew-1)*B_P(iELnew-1)/(B_P(iELnew-1)-B_P(iELnew))-Z_P(iELnew)*B_P(iELnew)/(B_P(iELnew-1)-B_P(iELnew));
            end
        end
        
        %Update some stuff
        CAPEallprev = CAPEall;
        CAPEst = sum(B_P(tropogridsize(2)+ mlgridsize(2)+1:iELnew).*dzst); %Actual calculated strato CAPE
        CAPEall = CAPEml + CAPEtp + CAPEst; %Actual calculated CAPE
        
        %New guess for stratosphere CAPE
        if abs(CAPEallprev - CAPEall)>=2
            CAPEstguess = (CAPEstguess+CAPEst)./2;
        else %if the CAPE iteration is stuck, force a pertubation
            CAPEstguess = CAPEst+(CAPEall-CAPE)*2;
        end
        
        %start trigger when calculated CAPE approach CAPE
        if abs(CAPEall-CAPE) < 1.5*CAPEblur
            ELguesstrigger = 1;
            ELguesstriggeraccum = ELguesstriggeraccum + ELguesstrigger;
        else
            ELguesstrigger = 0;
        end
        
        if abs(CAPEall-CAPE)<CAPEblur && CAPEendloopaccum<1
            CAPEendloop = 1;
            CAPEendloopaccum = CAPEendloopaccum + CAPEendloop;
        else
            CAPEendloop = 0;
             CAPEendloopaccum = 0;
        end
        
        ['still working on it... CAPE now at ',num2str(CAPEall)] %in case folks get anxious
    end
    
    if abs(ELnew-EL) > dzst
        ['EL unchievable - closest achievable EL chosen at ',num2str(ELnew),' m AGL']
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% below Zml %%%%%%%%%%%%%%%%%%%%%%%%
    
    CINtp = sum(B_P(mlgridsize(2)+1:iLFC).* dztp);
    Bml = g.*(Tv_P(mlgridsize(2)) - Tvml)/Tvml;
    
    syms a b
    eqns = [a*Zml*(exp(b*(Zml-LFC))-1) == Bml,...
        -(a*Zml^2)/2 + (a*exp(-b*LFC)*(1 + exp(b*Zml)*(-1 + b*Zml)))*b^(-2)== -CIN-CINtp];
    S =vpasolve(eqns,[a b],[-5e-5 -1e-3]);
    mla1 = double(S.a);
    mlb1 = double(S.b);
    
    if isempty(mla1) %if this formula does not work for our CIN, choose alternative
        ['alt buoyancy profile below Zml used']
        syms c d
        eqns2 = [c.*Zml.*(Zml-d) == Bml,...
            c*Zml^2*(2*Zml - 3*d)/6== -CIN-CINtp];
        S =vpasolve(eqns2,[c d]);
        mla2 = double(S.c);
        mlb2 = double(S.d);
        B_P(1:mlgridsize(2)) = mla2.*Z_P(1:mlgridsize(2)).*(Z_P(1:mlgridsize(2))-mlb2);
    else
        B_P(1:mlgridsize(2)) = mla1.*Z_P(1:mlgridsize(2)).*(exp(mlb1.*(Z_P(1:mlgridsize(2))-LFC))-1);
    end
    
    for i=1:mlgridsize(2)
        Tv_ml(i) = Tv_P(i).*g./(B_P(i)+g);
        T_ml(i)=  Tv_ml(i)./((1+rt_ml(i)./epsi)./(1+rt_ml(i)));
        
        Cp = Cpd*(1+rt_ml(i).*(Cpv./Cpd))./(1+rt_ml(i));
        R_ml(i) = Rd.*(1+rt_ml(i)./epsi)./(1+rt_ml(i));
        Theta_ml(i)=T_ml(i)*(P0./P_ml(i)).^(R_ml(i)./Cp);
        es_ml(i)=es(T_ml(i));
        H_ml(i)=rt_ml(i)*P_ml(i)./(es(T_ml(i))*(epsi+rt_ml(i)));
        
    end
end

%% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% compose output data matrix  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sdZ=transpose([Z_ml Z_tp Z_st]);
sdTheta=transpose([Theta_ml Theta_tp Theta_st]);
sdP=transpose([P_ml P_tp P_st]);
sdes=transpose([es_ml es_tp es_st]);
sdRH=transpose([H_ml H_tp H_st]);
sdT=transpose([T_ml T_tp T_st]);
sdqv=rt(sdP,sdes,sdRH);
sdTv=sdT.*(1+sdqv./epsi)./(1+sdqv);
sdTvP=transpose([Tv_P]);
sdB=g.*(sdTvP-sdTv)./sdTv;
sdNsign = (sdTheta(2:end)-sdTheta(1:end-1))./(sdZ(2:end)-sdZ(1:end-1));
CIN = -sum(sdB(1:find(sdB>1e-4,1))*25); %1e-4 m/s2 tolerance: B ~ 0 in the mixed layer when lapse_ml = g/Cpg

for umax=(30:5:30)  %umax at 6km
    for umid=(20:5:20)  %umax at 3km

        %%%%%%%%%%%%%%%%%%%%%%% wind speed and direction profile %%%%%%%%%%%%%%%%%%
        sdgridsize=alli;
        WD=zeros(sdgridsize,1);
        W=zeros(sdgridsize,1);
        sduv=makeloopy_fxcurve(ustrt,umid,umax,vstrt,vmax,sdP',sdZ',SRH);
        %sduv=makeloopy_allcross(umax,umid,sdP',sdZ');
        sdu=sduv(1,:)';
        sdv=sduv(2,:)';
        savnam = [ 'Sounding_Variables.mat'];
        save(savnam,'sdZ','sdTheta','sdP','sdes','sdRH','sdT','sdqv','sdu','sdTv','sdTvP','sdv','sdB','WD','W');
        savefig(['hodo_',num2str(CAPE),'_',num2str(umax),'_',num2str(umid),'_',num2str(SRH),'.fig'])
        print(['hodo_',num2str(CAPE),'_',num2str(umax),'_',num2str(umid),'_',num2str(SRH),'.png'],'-dpng')
         close
        %%%%%%%%%%%%%%%%%%%%% output data into txt sounding file  %%%%%%%%%%%%%%%%%

        %sfc pres (mb)    sfc theta (K)    sfc qv (g/kg)
        sfc=[sdP(1)./100; sdTheta(1); sdqv(1).*1000];
        %z (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)
        sddata=horzcat(sdZ, sdTheta, sdqv.*1000, sdu, sdv);

        %delete input_sounding

        fileID = fopen(['input_sounding_',num2str(umax),'_',num2str(umid)],'w');
        fprintf(fileID,'%12.4f %12.4f %12.4f\n',sfc);

        for i = 1:sdgridsize
            writedata=sddata(i,:);
            fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f\n',writedata);
        end


        fclose(fileID);


    end
end
close
tskew(sdTvP,sdZ,sdP./100,sdT-273.15,sdRH);
savefig(['skt_',num2str(CAPE),'.fig'])
print(['skt_',num2str(CAPE),'.png'],'-dpng')



function e=es(T)
e=6.112*exp(17.67*(T-273.15)./(T-273.15+243.15))*100;
end

function success=tskew(TvP,zz,pz,tz,rhz)

ez=es(tz+273.15);
qz=rt(pz*100,ez,rhz);
Tv=(tz+273.15).*(1+qz./0.622)./(1+qz)-273.15;
tvp=TvP-273.15;
chi=log(pz.*qz./(6.112.*(0.622+qz)));
tdz=243.5.*chi./(17.67-chi);
sz=size(pz);

 %% CAPE!


B=9.8.*(TvP-(Tv+273.15))./(Tv+273.15);

P1=1;
P2=sz(1);
P3=sz(1);

 for i=2:sz(1)

     if (i > P1) && ( B(i) >=0 ) && (B(i-1) <=0 )
         P2=i;
     end
     
     if (i > P2) && ( B(i) <=0 )  && (B(i-1) >=0 )
         P3=i;
         P3
     end
     
     
 end
 
CIN=trapz(zz(1:P2),B(1:P2));
CAPE = trapz(zz(P2:P3),B(P2:P3));


 %% Plot
 
figure(1)

%
p=[1025:-25:50];
pplot=transpose(p);
t0=[-40:2:40];
[ps1,ps2]=size(p);
ps=max(ps1,ps2);
[ts1,ts2]=size(t0);
ts=max(ts1,ts2);

tem=zeros(ts,ps);
thet=zeros(ts,ps);
tv=zeros(ts,ps);
q=zeros(ts,ps);
thetaea=zeros(ts,ps);

for i=1:ts
   for j=1:ps
      tem(i,j)=t0(i)+30.*log(0.001.*p(j));
      thet(i,j)=(273.15+tem(i,j)).*(1000./p(j)).^.287-273.15;
      tv(i,j)=thet(i,j)*(1+0.6*qz(j));
      e=6.112.*exp(17.67.*tem(i,j)./(243.5+tem(i,j)));
      q(i,j)=622.*e./(p(j)-e);
      thetaea(i,j)=thet(i,j).*exp(2.5.*q(i,j)./(tem(i,j)+273.15));
   end
end
p=transpose(p);
t0=transpose(t0);
temp=transpose(tem);
theta=transpose(thet);
thetae=transpose(thetaea);
tva=transpose(tv);
qs=transpose(sqrt(q));
h=contour(t0,pplot,temp,16,'k');
hold on
x0=10;
y0=10;
width=500;
height=600;
set(gcf,'position',[x0,y0,width,height])
%set(gca,'ytick',[1020:20:800])
set(gca,'yscale','log','ydir','reverse')
set(gca,'FontSize',12)
set(gca,'ytick',[100:100:1000])
%set(gca,'ytick',[800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020])
set(gca,'ygrid','on')
axis([-30 40 150 1000])
hold on
h=contour(t0,pplot,theta,24,'b');
h=contour(t0,pplot,tva,24,'c');
h=contour(t0,pplot,qs,24,'g');
h=contour(t0,pplot,thetae,24,'r');
%tsound=30.+43.5.*log(0.001.*p);
%tsoundm=tsound-30.*log(0.001.*p);
tzm=tz-30.*log(0.001.*pz);
tdzm=tdz-30.*log(0.001.*pz);
tvm=Tv-30.*log(0.001.*pz);
tvpplot=tvp-30.*log(0.001.*pz);
h2=plot(tzm,pz,'r',tdzm,pz,'g--',tvm,pz,'b',tvpplot,pz,':b');
hold off
legend(h2,'Temperature','Dew Point','Virtual Temperature','Air Parcel Virtual Temperature');
set(h2,'linewidth',2)
xlabel('T, T_d (deg C)','FontSize',15)
ylabel('Pressure (hPa)','FontSize',15)
title( ['CAPE : ', num2str(CAPE),' ; CIN:' num2str(CIN)] )

success='yes';
end

function res=lcl_resid(H,Tg,Pg,zsat,Rd,Cpd,Cpv,g,P0,epsi)
%residual: parcel saturation mixing ratio at zsat minus surface mixing ratio
%for surface RH = H, in a well-mixed (constant theta_v) layer
rtg = rt(Pg,es(Tg),H);
Rg = Rd.*(1+rtg./epsi)./(1+rtg);
Cpg = Cpd*(1+rtg.*(Cpv./Cpd))./(1+rtg);
Thetavg = Tg*((P0./Pg)^(Rg./Cpg))*(1+rtg/epsi)/(1+rtg);
kap = Rd/Cpd;
Psat = P0*((Pg/P0)^kap - g*zsat/(Cpd*Thetavg))^(1/kap);
res = rt(Psat,es(Tg-g./Cpg*zsat),1) - rtg;  %parcel follows its own dry adiabat g/Cpg
end

function r=rt(P,es,RH)

r=0.622.*es.*RH./(P-RH.*es);

end

function sduv=makeloopy_fxcurve(ustrt,umid,umax,vstrt,vmax,pinterp,zlev,SRH0)
% SRH = 260;
% ss = 1;
% pinterp=sdP';
% zlev=sdZ';
% % rot=-rot0;
%
% ustrt=0; %u at 0m
% vstrt=0; %v at 0m
% umax=31;
% vmax=7;

%%%%%%%%%%%%
b = vmax;
zlist = (50:100:20125);
numk = numel(zlist);
udep1=1000.0;
udep2=3000.0;
udep3=6000.0;

qu0 = zeros(numk,1);
qv0 = qu0;
u0 = zeros(numk,1);
v0 = u0;

%% test quarter circle

for k=1:numk

    if zlist(k) <= udep1
        angl = pi - pi./2 .* (zlist(k)-zlist(1))./(udep1- zlist(1));
        %make loop and stuck to y axis
        center = b; % how much we move the loop
        quloop0 = b.*cos(angl);
        qvloop0 = b.*sin(angl);
        quloop1 = quloop0+center;
        qvloop1 = qvloop0;
        %put the whole thing together
        utemp1 = quloop1;
        vtemp1 = qvloop1;

        %move the thing
        qu0(k) = utemp1+ustrt;
        qv0(k) = vtemp1+vstrt;
        if zlist(k+1) >  udep1
            qu00 = qu0(k);
            qv00 = qv0(k);
            qalluL = sum(pinterp(2:k+1).*qu0(1:k)');
            kl = k;
        end

    elseif zlist(k) > udep1 && zlist(k) <= udep3
        %make the straight part
        if zlist(k+1) >  udep2 && zlist(k) <= udep2
            km = k;
        end

        qu0(k) = qu00+(zlist(k)-udep1).*(umax-qu00)./(udep3-udep1) ;
        qv0(k) = qv00;

        if zlist(k+1) >  udep3
            kh = k;
            u22=qu0(k);
            v22=qv0(k);
        end

    elseif zlist(k) > udep1 && zlist(k) <= udep2

        u0(k) = qu00+(zlist(k)-udep1).*(umid-qu00)./(udep2-udep1);
        v0(k) = v00;


        if zlist(k+1) >  udep2
            u11=u0(k);
            v11=v0(k);
        end

    elseif zlist(k) > udep2 && zlist(k) <= udep3

        u0(k) = u11+(zlist(k)-udep2).*(umax-umid)./(udep3-udep2);
        v0(k) = v11;

        if zlist(k+1) >  udep3
            u22=u0(k);
            v22=v0(k);
        end

    else
        u0(k) = umax;
        v0(k) = v22;
    end

end

%
% figure(2)
% plot(qu0(1:find(zlist>=udep1,1,'first')),qv0(1:find(zlist>=udep1,1,'first')),'-ok','LineWidth',1)
% hold on
% plot(qu0(find(zlist>=udep1,1,'first'):find(zlist>udep2,1,'first')),qv0(find(zlist>=udep1,1,'first'):find(zlist>udep2,1,'first')),'-ob','LineWidth',1)
% plot(qu0(find(zlist>=udep2,1,'first'):find(zlist>udep3,1,'first')),qv0(find(zlist>=udep2,1,'first'):find(zlist>udep3,1,'first')),'-or','LineWidth',1)
% plot([-50,50],[0,0],'k')
% grid on
% xlim([-20 40])
% ylim([-20 20])



%% Bunkers attempt

%load('pinterp.mat')   %Standard atmosphere p/height relation
zlist1 = [0, zlist];
qu_raw = [0, qu0']; qv_raw = [0, qv0'];
%zlev = (0:10:20000);   %same as in standard atmosphere calculation

quinterp = interp1(zlist1,qu_raw,zlev);
qvinterp = interp1(zlist1,qv_raw,zlev);
qsduv(1,:)=quinterp;
qsduv(2,:)=qvinterp;

lt=find(zlev>500);
db=find(zlev>=5500);
dt=find(zlev>6000);
qlowu = mean(quinterp(1:lt(1)-1),"omitnan"); qlowv = mean(qvinterp(1:lt(1)-1),"omitnan");
qdeepu = mean(quinterp(db(1):dt(1)-1),"omitnan"); qdeepv = mean(qvinterp((db(1):dt(1)-1)),"omitnan");

%Pressure-weighted mean
qmeanu = mean(pinterp(1:dt(1)-1).*quinterp(1:dt(1)-1),"omitnan")./mean((pinterp(1:dt(1)-1)));
qmeanv = mean(pinterp(1:dt(1)-1).*qvinterp(1:dt(1)-1),"omitnan")./mean((pinterp(1:dt(1)-1)));


slope_orth = (qdeepv - qlowv)./(qdeepu - qlowu);
theta_1 = atand(slope_orth);

qdel_u = 7.5 .* sind(theta_1);
qdel_v = 7.5 .* cosd(theta_1);

%Right mover
RMu = qmeanu + qdel_u;
RMv = qmeanv - qdel_v;

%Left mover
LMu = qmeanu - qdel_u;
LMv = qmeanv + qdel_v;

%figure(2)
% plot([qlowu qdeepu],[qlowv qdeepv],'k')
% plot([qmeanu qmeanu],[qmeanv qmeanv],'x','MarkerSize',15)
% plot([RMu RMu],[RMv RMv],'*','MarkerSize',10)
% plot([LMu LMu],[LMv LMv],'*','MarkerSize',10)
% set(gcf,'position',[10,10,800,500])
%
qumove =RMu;
qvmove =RMv;

qSRFall=0;
for xx=1:find(zlist>=2000,1,'first')
    qSRFall = qSRFall + ((qu0(xx)-qumove).^2+(qv0(xx)-qvmove).^2).^0.5;
end

qSRF=qSRFall./find(zlist>=2000,1,'first');

qSRH=0;
for xx=1:find(zlist>=3000,1,'first')
    qSRH = qSRH + abs((qu0(xx+1)-qumove).*(qv0(xx)-qv0(xx+1))) + abs((qv0(xx)-qvmove).*(qu0(xx+1)-qu0(xx))) + abs((qv0(xx)-qv0(xx+1)).*(qu0(xx)-qu0(xx+1)));
end

%% ellipse shaped loopy hodograph

SRHg = 0;
br = 0.01;

while abs(SRHg - SRH0) >= 0.5
    abc = br.*b;
    for k=1:numk
        bnew = b - abc;
        alpha0 = pi - asin (-abc./bnew);
        anew = -b./cos(alpha0);
        %make loop and stuck to y axis
        center = b; % how much we move the loop

        if zlist(k) <= udep1

            vloop1 = qv0(k);
            uloop1 = cos(pi + asin((vloop1 - abc) ./ bnew)) .* anew + center;

            %put the whole thing together
            utemp1=uloop1;
            vtemp1=vloop1;

            %move the thing
            u0(k)=  utemp1;
            v0(k)=  vtemp1;
            if zlist(k+1) >  udep1
                u00=u0(k);
                v00=v0(k);
                alluL = sum(pinterp(2:k+1).*u0(1:k)');
                alluLdiff = (alluL - qalluL)./(km-kl-1);
            end


        elseif zlist(k) > udep1 && zlist(k) <= udep2

            u0(k) = qu00+(zlist(k)-udep1).*(umid-qu00)./(udep2-udep1);
            v0(k) = v00;


            if zlist(k+1) >  udep2
                u11=u0(k);
                v11=v0(k);
            end

        elseif zlist(k) > udep2 && zlist(k) <= udep3

            u0(k) = u11+(zlist(k)-udep2).*(umax-umid)./(udep3-udep2);
            v0(k) = v11;

            if zlist(k+1) >  udep3
                u22=u0(k);
                v22=v0(k);
            end

        else
            u0(k) = umax;
            v0(k) = v22;
        end

    end





    %% Bunkers attempt

    %load('pinterp.mat')   %Standard atmosphere p/height relation
    zlist1 = [0, zlist];
    u_raw = [0, u0']; v_raw = [0, v0'];
    %zlev = (0:10:20000);   %same as in standard atmosphere calculation

    uinterp = interp1(zlist1,u_raw,zlev);
    vinterp = interp1(zlist1,v_raw,zlev);
    sduv(1,:)=uinterp;
    sduv(2,:)=vinterp;

    lt=find(zlev>500);
    db=find(zlev>=5500);
    dt=find(zlev>6000);
    lowu = mean(uinterp(1:lt(1)-1),"omitnan"); lowv = mean(vinterp(1:lt(1)-1),"omitnan");
    deepu = mean(uinterp(db(1):dt(1)-1),"omitnan"); deepv = mean(vinterp((db(1):dt(1)-1)),"omitnan");

    %Pressure-weighted mean
    meanu = mean(pinterp(1:dt(1)-1).*uinterp(1:dt(1)-1),"omitnan")./mean((pinterp(1:dt(1)-1)));
    meanv = mean(pinterp(1:dt(1)-1).*vinterp(1:dt(1)-1),"omitnan")./mean((pinterp(1:dt(1)-1)));


    slope_orth = (deepv - lowv)./(deepu - lowu);
    theta_1 = atand(slope_orth);

    del_u = 7.5 .* sind(theta_1);
    del_v = 7.5 .* cosd(theta_1);

    %Right mover
    RMu = meanu + del_u;
    RMv = meanv - del_v;

    %Left mover
    LMu = meanu - del_u;
    LMv = meanv + del_v;

    umove =RMu;
    vmove =RMv;

    %%%%%%%%%%%%%%
    SRFall1=0;
    for xx=1:find(zlist>=1000,1,'first')
        SRFall1 = SRFall1 + ((u0(xx)-umove).^2+(v0(xx)-vmove).^2).^0.5;
    end
    SRF1=SRFall1./find(zlist>=1000,1,'first');

    SRFall2=0;
    for xx=1:find(zlist>=2000,1,'first')
        SRFall2 = SRFall2 + ((u0(xx)-umove).^2+(v0(xx)-vmove).^2).^0.5;
    end
    SRF2=SRFall2./find(zlist>=2000,1,'first');

    SRFall3=0;
    for xx=1:find(zlist>=3000,1,'first')
        SRFall3 = SRFall3 + ((u0(xx)-umove).^2+(v0(xx)-vmove).^2).^0.5;
    end
    SRF3=SRFall3./find(zlist>=3000,1,'first');
    %%%%%%%%%%%%%

    SRHg=0;
    for xx=1:find(zlist>=1000,1,'first')
        SRHg = SRHg + abs((u0(xx+1)-0).*(v0(xx)-v0(xx+1))) + abs((v0(xx)-0).*(u0(xx+1)-u0(xx))) + abs((v0(xx)-v0(xx+1)).*(u0(xx)-u0(xx+1)));

    end
    if SRHg < SRH0
        br = br + 0.0005;
    else
        br = br - 0.0005;
    end
end

SRH1 = 0;
for xx=1:find(zlist>=1000,1,'first')
    SRH1 = SRH1 + abs((u0(xx+1)-umove).*(v0(xx)-v0(xx+1))) + abs((v0(xx)-vmove).*(u0(xx+1)-u0(xx))) + abs((v0(xx)-v0(xx+1)).*(u0(xx)-u0(xx+1)));
end

SRH2 = 0;
for xx=1:find(zlist>=2000,1,'first')
    SRH2 = SRH2 + abs((u0(xx+1)-umove).*(v0(xx)-v0(xx+1))) + abs((v0(xx)-vmove).*(u0(xx+1)-u0(xx))) + abs((v0(xx)-v0(xx+1)).*(u0(xx)-u0(xx+1)));
end
%
% SRH1
% SRH2
SRH3 = 0;
for xx=1:find(zlist>=3000,1,'first')
    SRH3 = SRH3 + abs((u0(xx+1)-umove).*(v0(xx)-v0(xx+1))) + abs((v0(xx)-vmove).*(u0(xx+1)-u0(xx))) + abs((v0(xx)-v0(xx+1)).*(u0(xx)-u0(xx+1)));
end

figure(2)
plot(u0(1:find(zlist>=udep1,1,'first')),v0(1:find(zlist>=udep1,1,'first')),'-k','LineWidth',5)
hold on
plot(u0(find(zlist>=udep1,1,'first'):find(zlist>udep2,1,'first')),v0(find(zlist>=udep1,1,'first'):find(zlist>udep2,1,'first')),'-b','LineWidth',5)
plot(u0(find(zlist>=udep2,1,'first'):find(zlist>udep3,1,'first')),v0(find(zlist>=udep2,1,'first'):find(zlist>udep3,1,'first')),'-r','LineWidth',5)
plot([-50,50],[0,0],'k')
plot([0,0],[-50,50],'k')
grid on
xlim([-11 35])
ylim([-5 15])

figure(2)
plot([lowu deepu],[lowv deepv],'k')
plot([meanu meanu],[meanv meanv],'x','MarkerSize',15)
plot([RMu RMu],[RMv RMv],'*','MarkerSize',10)
plot([LMu LMu],[LMv LMv],'*','MarkerSize',10)
set(gcf,'position',[10,10,800,400])

text(-10,12.5,{['SRH01, 02, 03 = ',num2str(SRH1,'%.2f'),', ',num2str(SRH2,'%.2f'),', ',num2str(SRH3,'%.2f'),' m^2 s^-^2']...
    ['SRF01, 02, 03 = ',num2str(SRF1,'%.2f'),', ',num2str(SRF2,'%.2f'),', ',num2str(SRF3,'%.2f'),' m s^-^1']...
    ['Bunkers u,v = ',num2str(RMu,'%.2f'),',',num2str(RMv,'%.2f')]},'FontSize',10)

end