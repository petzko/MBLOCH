function plotQCL(Psi,E,Vx,x,wavefunc_flag,modulisqr_flag,modulisqrtrunc_flag,heading)

nx = length(Psi(:,1));
nlevel = length(E);
xA = x*1e10;

PsiPsi = zeros(nx,nlevel);

for ii = 1 : nlevel
    for kk = 1 : nx
        PsiPsi(kk,ii) = abs(Psi(kk,ii)) * abs(Psi(kk,ii));
    end
end

legendInfo= cell(nlevel+1,1); 
legendInfo{1} = 'V_c';


for ii = 1 : nlevel
    PsiPlot(:,ii) = Psi(:,ii)/(5*ii*max(Psi(:,ii))) + E(ii);
    legendInfo{ii+1} = ['|\Psi_{' num2str(ii) '}|^2'];
end


if wavefunc_flag == 1
    figure;
    plot(xA,Vx,'black','LineWidth',2)
    plot(xA,PsiPlot,'LineWidth',2);
    legend(legendInfo);
    set(gca,'fontname','Times New Roman','fontsize',20);
    set(gca,'LineWidth',1.5,'TickLength',[0.03 0.025]); 
    xlabel('Position (Angstrom)','fontname','Times New Roman','fontsize',20);
    ylabel('Energy (eV)','fontname','Times New Roman','fontsize',20);
    hold on;

    set(gcf, 'WindowStyle', 'normal');
    set(gca, 'Unit', 'inches');
    set(gca, 'Position', [0.9 0.75 4.2 3.3]);   
    %title('\Psi','FontSize',16);
    
end

if modulisqr_flag == 1
    for ii = 1 : nlevel
        PsiPsiPlot(:,ii) = PsiPsi(:,ii)/(20*max(PsiPsi(:,ii))) + E(ii);
    end

    figure;
    plot(xA,Vx,'black','linewidth',2);
    hold on;
    plot(xA,PsiPsiPlot,'linewidth',2);
    legend(legendInfo);
%     set(gca,'fontname','Times New Roman','fontsize',20);
%     set(gca,'LineWidth',1.5,'TickLength',[0.03 0.025]);
    xlabel('Position (Angstrom)','fontname','Times New Roman','fontsize',20);
    ylabel('Energy (eV)','fontname','Times New Roman','fontsize',20);
  
    %hold on;
%     axis 'tight';
%     set(gcf, 'WindowStyle', 'normal');
%     set(gca, 'Unit', 'inches');
%     set(gca, 'Position', [0.9 0.75 4.2 3.3]);   
    %title('\Psi* \Psi','FontSize',16);
    
end

if modulisqrtrunc_flag
    figure;
    plot(xA,Vx,'black','LineWidth',2)
    set(gca,'fontname','Times New Roman','fontsize',20);
    set(gca,'LineWidth',1.5,'TickLength',[0.03 0.025]);
    xlabel('Position (Angstrom)','fontname','Times New Roman','fontsize',20);
    ylabel('Energy (eV)','fontname','Times New Roman','fontsize',20);
    hold on;
set(gcf, 'WindowStyle', 'normal');
set(gca, 'Unit', 'inches');
set(gca, 'Position', [0.9 0.75 4.2 3.3]);   
    colors = get(gca,'ColorOrder');
    color = 1;
    levelindex = 1;
    for ii = 1 : nlevel
    
        if ii < 8
            color = ii;
        else
            if ii < 15
                color = ii - 7;
            else
                if ii < 22
                    color = ii - 14;
                else
                    if ii < 29
                        color = ii - 21;
                    else
                        color = ii - 28;
                    end
                
                end
            end
        end

    
        jj = 0;
        sindex_flag = 1;
        findex_flag = 0;
    
        for kk = 1 : nx
            if PsiPsiPlot(kk,ii) > (E(ii) + 5e-3)
        
                jj = jj + 1;
    
                if jj == 1
                    findex_flag = 1;
                end
            
                if sindex_flag == 1
                    sindex = kk;
                    sindex_flag = 0;
                end

            else
                if findex_flag == 1
                    findex = kk - 1;
                    findex_flag = 0;
                end
            end   
            jj = 0;
        end
        selectcolor = colors(color,:);
        plot(xA(sindex:findex),PsiPsiPlot(sindex:findex,ii),'color',selectcolor,'LineWidth',2);
        selectcolor = colors(color);
        hold on;
    end
    %title('\Psi* \Psi','FontSize',16);
end
      title(heading);