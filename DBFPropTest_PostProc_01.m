clear all;
close all;

% Establish working filename(s) and directory...
rootpath = pwd;
selregsavfname = [rootpath '\selregsave.mat'];

% Set test constants...
Tamb = (66+460);
RH = 74;
pamb = 28.81*(2116/29.92);

% Determine ambient air density with RH correction...
es = 1.7526e11*exp(-5315.56*1.8/Tamb);
rhoambSI = (0.0034847*1.8/Tamb)*(pamb*(101325/2116)-0.003796*RH*es);
rhoamb = rhoambSI*(0.002377/1.225);

% Get run parameters and find test section velocity...
matrixfname = 'DBF_PropTest_01.txt';
FID = fopen([rootpath '\' matrixfname]);
mfiledata = cell2mat(textscan(FID,'%f%f%f%f%f%f','HeaderLines',1,'Delimiter','\t'));
fclose(FID);
mfiledata(:,end+1) = (2*5.2*mfiledata(:,6)/rhoamb).^0.5;

% Get recorded data from RCBenchmark DAQ and perform calcs...
listing = dir(rootpath);
fnamect = 0;
dfnamect = 0;
for i = 1:1:length(listing)
    if ~listing(i).isdir
        fnamect = fnamect + 1;
        fnameind(fnamect) = i;
        if listing(i).name(end-3:end) == '.csv'
            dfnamect = dfnamect + 1;
            dfnameind(dfnamect) = i;
        end
    end
end

% Check for previous extraction of stable RPM region...
if exist(selregsavfname, 'file') == 2
    ldsw = menu('Use previously-extracted regions of stable RPM?','Yes','No');
    if ldsw == 1
        load(selregsavfname,'bracket');
    end
else
    ldsw = 2;
end

for i = 1:1:length(dfnameind)
    % Setup run parameters...
    DEng_ft = mfiledata(i,2)/12;
    VTS_fps = mfiledata(i,7);
    
    % Get collected run data...
    rundata = importdata([rootpath '\' listing(dfnameind(i)).name]);
    rdata{i} = rundata.data;
    t_sec{i} = rdata{i}(:,1);
    ESCsig{i} = rdata{i}(:,2);
    T_lbf{i} = rdata{i}(:,10);
    V_VDC{i} = rdata{i}(:,11);
    i_A{i} = rdata{i}(:,12);
    Om_RPM{i} = rdata{i}(:,13);
    vib_g{i} = rdata{i}(:,20);
    
    % Perform calcs...
    P_ftlbfpsec{i} = (550/746)*V_VDC{i}.*i_A{i};
    n_rps{i} = Om_RPM{i}/60;
    CT_raw{i} = T_lbf{i}./(rhoamb*DEng_ft^4*n_rps{i}.^2);
    CP_raw{i} = P_ftlbfpsec{i}./(rhoamb*DEng_ft^5*n_rps{i}.^3);
    if VTS_fps==0
        J_raw{i} = NaN*ones(size(t_sec{i}));
        eta_raw{i} = NaN*ones(size(t_sec{i}));
    else
        J_raw{i} = VTS_fps./(n_rps{i}*DEng_ft);
        eta_raw{i} = CT_raw{i}.*J_raw{i}./CP_raw{i};
    end
    
    % Allow user to manually select region of stable RPM...
    if ldsw ~= 1
        plot(n_rps{i},'.');
        disp('Select region of stable RPM for current test');
        bracket{i} = ginput(2);
    end
    selind = round(bracket{i}(:,1));
    CT{i} = CT_raw{i}(min(selind):max(selind));
    CP{i} = CP_raw{i}(min(selind):max(selind));
    J{i} = J_raw{i}(min(selind):max(selind));
    Om{i} = Om_RPM{i}(min(selind):max(selind));
    eta{i} = eta_raw{i}(min(selind):max(selind));
    close all;
    clear rundata
end

% Plot results...
figno = 0;
for i = 1:1:length(dfnameind)
    titname = sprintf('%i x %i Propeller, V_{TS} = %0.2f ft/sec',....
        mfiledata(i,2),mfiledata(i,3),mfiledata(i,7));
    if mfiledata(i,7)==0
        Xdata = Om{i};
        Xlab = '\Omega (RPM)';
    else
        Xdata = J{i};
        Xlab = 'J';
        Ydata3 = eta{i};
        figno = figno + 1;
        figure(figno);
        plot(Xdata,Ydata3,'gd');
        xlabel(Xlab); ylabel('\eta');
        grid on; grid minor;
        title(titname);
        ylim([0 1]);
    end
    Ydata1 = CT{i};
    Ydata2 = CP{i};
    legname = {'C_P', 'C_T'};
    figno = figno + 1;
    figure(figno);
    plot(Xdata,Ydata1,'ks',Xdata,Ydata2,'ro');
    xlabel(Xlab);
    grid on; grid minor;
    legend(legname,'Location','Best');
    title(titname);
end

% Allow the user to save the selected regions of stable RPM...
if ldsw ~= 1
    svsw = menu('Save selection of stable RPM regions?','Yes','No');
    if svsw == 1
        save([rootpath '\selregsave.mat'],'bracket');
    end
end

% Allow the user the option to save the plots...
prsw = menu('Save figures?','Yes','No');
if prsw == 1
    for i = 1:1:figno
        figure(i);
        print([rootpath '\' sprintf('fig_%02d.jpg',i)],'-djpeg','-r300');
    end
end