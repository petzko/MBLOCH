function output = SWIFTS_TD(data,noise_range,plot_ranges,filter_time)
% data: structure from Plot_SWIFTS
% noise_range: 1x2 vector specifying the noise range in THz (used to find
%              std deviation of measurement)
% plot_range:  Mx2 matrix of ranges to plot
% filter_time: in ps (Hanning window)

NUM_INSTANCES = 1000;
SIGMAS = [-2,-1,0,1,2];

ts = data.ts-mean(data.ts);
fs = data.pk_fs;
cs = data.pk_cs;

gis = find(fs >= noise_range(1) & fs <= noise_range(2));
ns = cs(gis,:); ns = ns-repmat(mean(ns,1),size(ns,1),1);
ns_std = sqrt(0.5*mean(real(ns).^2 + imag(ns).^2,1))
% the noise in each quadrature is a zero-mean gaussian with std
% deviation ns_std
ns_insts = randn(size(cs,1),size(cs,2),NUM_INSTANCES)+i*randn(size(cs,1),size(cs,2),NUM_INSTANCES);
ns_insts = ns_insts.*repmat(ns_std,size(cs,1),1,NUM_INSTANCES);

nsy_cs  = repmat(cs,1,1,NUM_INSTANCES)-ns_insts;
Sps = (data.C*nsy_cs(:,1,:) - nsy_cs(:,2,:))/(i*imag(data.C));
Sms = (conj(data.C)*nsy_cs(:,1,:) - nsy_cs(:,2,:))/(-i*imag(data.C));
Sns = nsy_cs(:,3,:);
Sps=squeeze(Sps); Sms=squeeze(Sms); Sns=squeeze(Sns);

for wl=5

    fp  = fs(1:end-1);
    cp  = abs(Sps(1:end-1,:));
    cpp = sqrt(abs(Sns(1:end-1,:)).*abs(Sns(2:end,:)));


    for i1=1:NUM_INSTANCES
%         cpp(:,i1) = cpp(:,i1) * ((cpp(:,i1)'*cpp(:,i1))\(cpp(:,i1)'*cp(:,i1)));
        [~,ml]=max(cpp(:,i1));
        cpp(:,i1) = cpp(:,i1)*cp(ml,i1)/cpp(ml,i1);
    end

    w = hanning(wl); w=w/sum(w);
    M = length(fp)-wl+1;
    ls = repmat([1:wl]',1,M)+repmat([0:M-1],wl,1);                      % index matrix for degree of coherence
    ls2= repmat(ls,1,NUM_INSTANCES)+kron([0:NUM_INSTANCES-1]*length(fp),ones(size(ls)));
    w2 = repmat(w,1,M*NUM_INSTANCES);
    g = sum(cpp(ls2).*w2.*cp(ls2),1)./sum(cpp(ls2).*w2.*cpp(ls2),1);
    g = reshape(g,[M,NUM_INSTANCES]);
    
    fc = conv2(fp,w,'valid');
    g = sort(g,2);
    cis = round(1/2*(1+erf(SIGMAS/sqrt(2))) * (NUM_INSTANCES-1) + 1);
    
%     dfigure; set(gcf,'Position',[191 566 1470 222]); movegui('center');
    subplot(3,1,3);
    co = get(gca,'ColorOrder');
    plot(fc,g(:,cis(3)),'Color',co(1,:));
    hold all;
    pobj = Plot_Range(fc,g(:,cis(1)),g(:,cis(5)),co(1,:));
    ylim([0,1.5]);
end

output = data;
end

function pobj = Plot_Range(x,ybtm,ytop,color)
    xs = [x;flipud(x)];
    ys = [ybtm;flipud(ytop)];
    pobj = fill(xs,ys,color);
    set(pobj,'EdgeColor','none','FaceAlpha',0.2);
end