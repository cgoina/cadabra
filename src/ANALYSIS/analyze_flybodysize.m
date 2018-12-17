%% Analyze_flybodysize
% Copyright (C) 2008 Heiko Dankert, California Institute of Technology

% This file is part of ANALYSIS
% and the "Caltech Automated Drosophila Aggression-Courtship 
% Behavioral Repertoire Analysis (CADABRA)".

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.

% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% * Implementation by Heiko Dankert
%
%% Analyze and compare body length and area of flies

function FID = analyze_flybodysize(obj1,obj2,lunges,nmovs,gen_ident,gen_ident1,params,FID)

if sum(obj1{1}.am) > 0 && params.plots.stat.body,
    % FLY BODY AREA AND LENGTH
    ngens = length(lunges);
    utest = utestpairs(gen_ident1);
    x = 1:ngens; ya.mea = zeros(ngens,2); ya.sem = ya.mea; ya.min = ya.mea; ya.max = ya.mea; yl = ya;
    Xa = []; Ga = []; Xl = []; Gl = []; Xas = []; Gas = []; Xls = []; Gls = [];
    data1.area = cell(ngens,1); data1.length = data1.area; 
    data2.area = data1.area; data2.length = data1.area;
    for igen=1:ngens,
        tmp1 = []; tmp2 = []; tmp1l = []; tmp2l = [];
        for j=1:nmovs(igen),
            % Distinguish between winner (fly with more lunges) and loser?
            % Collect fly body areas (am) and lengths (lm) for fly 1 and
            % fly 2 or for winner and loser
            if (lunges{igen}.obj1.number(j) >= lunges{igen}.obj2.number(j)) || ~params.winlos || params.courtship,
                tmp1 = [tmp1 obj1{igen}.am(j)]; tmp2 = [tmp2 obj2{igen}.am(j)];
                tmp1l = [tmp1l obj1{igen}.lm(j)]; tmp2l = [tmp2l obj2{igen}.lm(j)];
            else
                tmp1 = [tmp1 obj2{igen}.am(j)]; tmp2 = [tmp2 obj1{igen}.am(j)];
                tmp1l = [tmp1l obj2{igen}.lm(j)]; tmp2l = [tmp2l obj1{igen}.lm(j)];
            end
        end
        
        % Cell arays for text output
        data1.area{igen} = tmp1; data1.length{igen} = tmp1l;
        data2.area{igen} = tmp2; data2.length{igen} = tmp2l;
        
        % Prepare fly body area and length vectors (X) and their
        % corresponding label vector (G) for box plots
        Xa = [Xa tmp1 tmp2];
        Ga = [Ga ones(1,length(tmp1))*(2*igen-1) ones(1,length(tmp2))*2*igen];
        Xl = [Xl tmp1l tmp2l];
        Gl = [Gl ones(1,length(tmp1l))*(2*igen-1) ones(1,length(tmp2l))*2*igen];

        % Compute mean (mea), standard error of the mean (sem), 
        % minimum (min) and maximum (max) of fly body areas (ya) and
        % lengths (yl)
        ya.mea(igen,:) = [mean(tmp1) mean(tmp2)];
        ya.sem(igen,:) = [std(tmp1) std(tmp2)]./sqrt(nmovs(igen));
        ya.min(igen,:) = [min(tmp1) min(tmp2)]; ya.max(igen,:) = [max(tmp1) max(tmp2)];
        yl.mea(igen,:) = [mean(tmp1l) mean(tmp2l)];
        yl.sem(igen,:) = [std(tmp1l) std(tmp2l)]./sqrt(nmovs(igen));
        yl.min(igen,:) = [min(tmp1l) min(tmp2l)]; yl.max(igen,:) = [max(tmp1l) max(tmp2l)];

        % Compute mean fly body area and length of fly 1+2 or winner+loser
        data = (tmp1+tmp2)/2; datal = (tmp1l+tmp2l)/2; 
        Xas = [Xas data]; Gas = [Gas ones(1,length(data))*igen];
        Xls = [Xls datal]; Gls = [Gls ones(1,length(datal))*igen];

        yas.mea(igen,:) = mean(data);
        yas.sem(igen,:) = std(data)./sqrt(nmovs(igen));
        yas.min(igen,:) = min(data); yas.max(igen,:) = max(data);
        yls.mea(igen,:) = mean(datal);
        yls.sem(igen,:) = std(datal)./sqrt(nmovs(igen));
        yls.min(igen,:) = min(datal); yls.max(igen,:) = max(datal);
    end
    
    % Prepare genotype pairs to be statistically tested
    utest = utestpairs(gen_ident1);

    % Bar plots or box plots and text output
    if numel(find(params.plots.stat.body_vec == 1)),
        % Fly length
        % Text output
        if params.bool_xls,
            textoutput([data1.length , data2.length],gen_ident,'length [mm]','Fly body',params);
        end
        % Plots
        ytit = 'fly body length [mm]'; titl = 'Fly Body Length';
        FID = plot_bar_msem(yl,utest,Xl,Gl,3,ytit,titl,FID,gen_ident,params);
        FID = plot_bar_msem(yls,utest,Xls,Gls,3,ytit,titl,FID,gen_ident,params);
    end
    if numel(find(params.plots.stat.body_vec == 2)),
        % Fly area
        % Text output
        if params.bool_xls,
            textoutput([data1.area , data2.area],gen_ident,'area [mm2]','Fly body',params);
        end
        % Plots
        ytit = 'fly body area [mm^2]'; titl = 'Fly Body Area';
        FID = plot_bar_msem(ya,utest,Xa,Ga,2,ytit,titl,FID,gen_ident,params);
        FID = plot_bar_msem(yas,utest,Xas,Gas,2,ytit,titl,FID,gen_ident,params);
    end

    
    % ANALYZE FLY BODY SIZE DIFFERENCE
    if ~params.oneobj,
        x = 1:ngens; ya.mea = zeros(ngens,1); ya.sem = ya.mea; ya.min = ya.mea; ya.max = ya.mea; yl = ya;
        Xa = []; Ga = []; Xl = []; Gl = []; clear('data');
        data.area = cell(ngens,1); data.length = data.area; 
        for igen=1:ngens,
            tmp1 = []; tmp2 = [];
            for j=1:nmovs(igen),
                % Relate body area and length: fly 1 over fly 2 or winner over loser
                if (lunges{igen}.obj1.number(j) >= lunges{igen}.obj2.number(j)) || ~params.winlos || params.courtship,
                    tmp1 = [tmp1 obj1{igen}.am(j) / obj2{igen}.am(j)];
                    tmp2 = [tmp2 obj1{igen}.lm(j) / obj2{igen}.lm(j)];
                else
                    tmp1 = [tmp1 obj2{igen}.am(j) / obj1{igen}.am(j)];
                    tmp2 = [tmp2 obj2{igen}.lm(j) / obj1{igen}.lm(j)];
                end
            end
            % Convert into percentage; prepare fly body area and length vectors (X) and their
            % corresponding label vector (G) for box plots
            tmp1 = (1-tmp1)*100; tmp2 = (1-tmp2)*100;
            Xa = [Xa tmp1]; Ga = [Ga ones(1,length(tmp1))*igen];
            Xl = [Xl tmp2]; Gl = [Gl ones(1,length(tmp2))*igen];

            ya.mea(igen) = mean(tmp1);
            ya.sem(igen) = std(tmp1)./sqrt(nmovs(igen));
            ya.min(igen) = min(tmp1); ya.max(igen) = max(tmp1);
            yl.mea(igen) = mean(tmp2);
            yl.sem(igen) = std(tmp2)./sqrt(nmovs(igen));
            yl.min(igen) = min(tmp2); yl.max(igen) = max(tmp2);                        
        
            % Cell arays for text output
            data.area{igen} = tmp1; data.length{igen} = tmp2;
        end

        % Prepare genotype pairs to be statistically tested
        utest = utestpairs(gen_ident1);

        % Bar plots or box plots
        if numel(find(params.plots.stat.body_vec == 3)),        
            % Text output
            if params.bool_xls,
                textoutput(data.length,gen_ident,'length diff [percent]','Fly body',params);
            end
            % Plots
            if params.winlos, ytit = 'winner <-> loser fly length [%]'; else ytit = 'fly1 - fly2 body length [%]'; end
            titl = 'Fly Body Length Difference';
            FID = plot_bar_msem(yl,utest,Xl,Gl,[-50 50],ytit,titl,FID,gen_ident,params);
        end
        if numel(find(params.plots.stat.body_vec == 4)),
            % Text output
            if params.bool_xls,
                textoutput(data.area,gen_ident,'area diff [percent]','Fly body',params);
            end
            % Plots
            if params.winlos, ytit = 'winner <-> loser fly area difference [%]'; else ytit = 'fly1 - fly2 body area [%]'; end
            titl = 'Fly Body Area Difference';
            FID = plot_bar_msem(ya,utest,Xa,Ga,[-50 50],ytit,titl,FID,gen_ident,params);
        end
    end
end

end