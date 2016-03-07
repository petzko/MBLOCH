dir = pwd;


%9
% lvls = [6 8 5 1]

% for TB approximation
%10 10.1 10.2 10.3 10.4
lvls =  [10 8 9 7 6 5 3]

%for decoupled we have @ 10.2
% lvls =  [11 9 10 8 7 6 4]


INJ_1 = 1; INJ_2=2; ULL  =3 ; LLL_1 = 4; LLL_2 = 5; DEPOP_1 = 6; DEPOP_2 =7;


% 11
%lvls = [9 8 6 4];
% 12
%lvls = [9 7 6 4];
%13
% lvls = [8 6 5 3];


Tmtx = zeros(length(lvls));
Wmtx = Tmtx;


%get the energies..
Ens = load('E_BOUND');
[sortedEns, order]  = sort(Ens(1:end,2));
sortedNrs = Ens(order,1);
Energy_lvls = sortedNrs(lvls);
Energies = sortedEns(lvls);

Energy_lvls

for i=1:length(Energy_lvls)
    for j = 1:length(Energy_lvls)
        if j~=i
         
            global sf;
            scan(dir,Energy_lvls(i), Energy_lvls(j), 1:5,[],3);
            Tmtx(i,j) = 1/sf;
            Wmtx(i,j) = sf;
       
        end
    end
end
%%
inj = [INJ_1 INJ_2]; ull = [ULL]; lll = [LLL_1,LLL_2]; depop = [DEPOP_1,DEPOP_2];

aggregate_levels = {inj,ull,lll,depop};

%inj_idx = 1; ull_idx = 2; lll_idx = 3; depop_idx =  4;

idx = [1:4];
w4 = zeros(4);

for i = idx
    for j = idx;
        if i~=j
            for k = aggregate_levels{i}
                for m = aggregate_levels{j}
                    w4(i,j) = w4(i,j) + Wmtx(k,m);
                end
            end
            
        end
        
    end
    
end


%%
display(['E1p = ' num2str(Energies(INJ_1))]);
display(['E0p = ' num2str(Energies(INJ_2))]);
display(['E3 = ' num2str(Energies(ULL))]);
display(['E2 = ' num2str(Energies(LLL_1))]);
display(['E2p = ' num2str(Energies(LLL_2))]);
display(['E1 = ' num2str(Energies(DEPOP_1))]);
display(['E0 = ' num2str(Energies(DEPOP_2))]);
display(' ');
display(['Energies = ' num2str(Energies') ]);
display(' ');

%%
display(' ');
display(['W_inj1 = ' num2str(Wmtx(INJ_1,:)/1E12)]);
display(['W_inj2 = ' num2str(Wmtx(INJ_2,:)/1E12)]);
display(['W_ull = ' num2str(Wmtx(ULL,:)/1E12)]);
display(['W_lll1 = ' num2str(Wmtx(LLL_1,:)/1E12)]);
display(['W_lll2 = ' num2str(Wmtx(LLL_2,:)/1E12)]);
display(['W_depop1 = ' num2str(Wmtx(DEPOP_1,:)/1E12)]);
display(['W_depop2 = ' num2str(Wmtx(DEPOP_2,:)/1E12)]);
display(' ');
display('aggregates!!!');
display(' ');
display(['W_inj = ' num2str(w4(1,:)/1E12)]);
display(['W_ull = ' num2str(w4(2,:)/1E12)]);
display(['W_lll = ' num2str(w4(3,:)/1E12)]);
display(['W_depop = ' num2str(w4(4,:)/1E12)]);


display(['Wmtx = ' num2str(reshape(Wmtx',1,size(Wmtx,1)*size(Wmtx,2))/1E12)]);




