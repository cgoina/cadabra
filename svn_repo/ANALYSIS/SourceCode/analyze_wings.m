%% Analyze_wings
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
%% Plot wing statistics and analyze fly orientation in case of a wing
%% extention of the opponent including time series of fly distances, wing
%% extension duration

function FID = analyze_wings(lunges,wings,copu,dists,dts,nframes,dxy,path,sorted_names,gen_ident,genotypes,params,chamber,FID)


% PLOT HEATMAPS AND RASTER PLOTS
FID = plot_wing_heatm_raster(lunges,wings,dts,nframes,dxy,gen_ident,params,chamber,FID);

% PLOT WING EXTENSION STATISTICS
labelfntsize = params.axisfontsize;
% Occurences or time spent?
if params.tim,
    maxy = params.range.stat.tim.wing;
else
    maxy = params.range.stat.occ.wing;
end
addname = ['analysis' params.slash];

% Summed over both flies?
bool_sumobj = 1;
ngens = max(genotypes); %k = factor(ngens);
utest = utestpairs(gen_ident);

if params.plots.stat.leftwing && numel(find(params.plots.stat.wing_vec == 1)),
    titl = 'Left Wing Extensions'; ylab = 'left wing extended';
    FID = plot_loser_winner_stat(wings.ext.l,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
end
if params.plots.stat.rightwing && numel(find(params.plots.stat.wing_vec == 1)),
    titl = 'Right Wing Extensions'; ylab = 'right wing extended';
    FID = plot_loser_winner_stat(wings.ext.r,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest);
end
if params.plots.stat.onewing && ~params.plots.movieclips && numel(find(params.plots.stat.wing_vec == 1)),
    titl = 'One Wing Extensions'; ylab = 'one wing extended';
    FID = plot_loser_winner_stat(wings.ext.l,FID,maxy,gen_ident,params,ylab,dts,titl,bool_sumobj,lunges,utest,wings.ext.r);
end

% RELATIVE FLY ORIENTATION IN CASE OF WING-EXTENSION
if params.plots.heatmaps.wingrel && min(params.plots.heatmaps.wingrel_vec) < 4,
    ngens = length(lunges);
    extfli = 'ext'; % wing extensions (ext) or wing flicks (fli)
    maxx = 10; % chamber.width, range of distances of interest
    max_int_manscale = 400; % max. occurence in histogram
    add = 0; % add. n frames before/after each sequence
    tstep = 0; % >0 only for heatmap-timeseries
    % phi0.min = 40; phi0.max = 60; % borders for wing angle (rel. to major body-axis)
    phi0.min = 70; phi0.max = 90; % borders for wing angle (rel. to major body-axis)
    lrb_vec = ['l' 'r'];
    shift_vec = [-1 0 1]; % extract period before, during, or after wing extension
    min_vortho = 0.1; %[mm/frame]
    min_vparal = 0.; %[mm/frame]
    min_vphi = 0. * dts{1}(1); %[deg/frame]

    xaxis = -maxx:5:maxx;
    xtick_lab = cell(1,length(xaxis)); for i=1:length(xaxis), xtick_lab{i} = num2str(-maxx+(i-1)*5); end
    rotdir1 = zeros(numel(shift_vec),numel(lrb_vec),ngens,2); rotdir2 = rotdir1;
    imgr = cell(ngens,numel(shift_vec)); imgl = imgr;
    for ilrb=1:2, % both left and right wing
        lrb = lrb_vec(ilrb);
        for ishift=2:2, % period during wing extension
            shift = shift_vec(ishift);
            % EXTRACT RELATIV ORIENTATION AND POSITION OF ONE FLY TOWARDS 
            % THE OTHER FLY, WHILE THE OTHER FLY IS EXTENDING A WING
            fn = [path addname num2str(ngens) '_gens_wing_' extfli '_' lrb '_add_' ...
                num2str(add) '_minmaxphi_' num2str(phi0.min) '_' num2str(phi0.max) '_frames'];
            if shift<0,
                fn = [fn '_before_ob_orient'];
            elseif shift>0
                fn = [fn '_after_ob_orient'];
            else
                fn = [fn '_ob_orient'];
            end
            ngens = max(genotypes);
            o1 = cell(1,ngens); o2 = o1;
            fid = fopen([fn '.mat'],'r');
            if (fid < 0) || params.analyze_new,
                h1 = waitbar(0,'Analyzing Data');
                for igen=1:ngens,
                    waitbar(igen/ngens,h1,gen_ident{igen});
                    dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
                    [ob1,ob2] = extract_wing_orient(igen,wings,extfli,lrb,round(add/dt),shift,phi0,genotypes,sorted_names,copu);
                    o1{igen} = ob1; o2{igen} = ob2;
                end
                close(h1);
                save(fn,'o1','o2');
            else
                load(fn);
            end

            % RELATIVE FLY POSITION HEATMAPS FOR WING EXTENSION PHASES
            minx = -maxx; miny = -maxx; maxy = maxx; dv = maxx/15;
            xv = minx:dv:maxx; yv = miny:dv:maxy;

            ngens = max(genotypes); %k = factor(ngens);
            max_int = 0;
            phi_velo1 = cell(ngens,2); phi_velo2 = phi_velo1;
            for igen=1:ngens,
                img = [];
                if ~params.courtship || params.oneobj,
                    [img,phi_velo1,rotdir] = analyze_wing_mov(o1,phi_velo1,img,igen,add,tstep,min_vortho,min_vparal,min_vphi,xv,yv,chamber,params.oneobj);
                    if numel(rotdir), rotdir1(ishift,ilrb,igen,:) = rotdir; end
                end
                if ~params.courtship || ~params.oneobj,
                    [img,phi_velo2,rotdir] = analyze_wing_mov(o2,phi_velo2,img,igen,add,tstep,min_vortho,min_vparal,min_vphi,xv,yv,chamber,params.oneobj);
                    if numel(rotdir), rotdir2(ishift,ilrb,igen,:) = rotdir; end
                end
                if strcmp(lrb,'r'), imgr{igen,ishift} = img; end
                if strcmp(lrb,'l'), imgl{igen,ishift} = img; end
                % Prepare circling direction information
                sumdir = sum(abs(rotdir1(ishift,ilrb,igen,:))+abs(rotdir2(ishift,ilrb,igen,:)));
                rotdir1(ishift,ilrb,igen,:) = rotdir1(ishift,ilrb,igen,:) / sumdir;
                rotdir2(ishift,ilrb,igen,:) = rotdir2(ishift,ilrb,igen,:) / sumdir;
                if max(max(img)) > max_int, max_int = ceil(max(max(img))); end
            end
        end
    end
    
    if ~params.autoscale || max_int <= 0 || isnan(max_int), 
        max_int = max_int_manscale; 
    end
    
    if numel(find(params.plots.heatmaps.wingrel_vec == 1)),
        % PLOT HEATMAPS OF RELATIVE FLY POSTION WHILE EXTENDING A WING
        ngens = max(genotypes);
        for ilrb=1:2,
            lrb = lrb_vec(ilrb);
            for ishift=2:2,
                shift = shift_vec(ishift);
                if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
                figure(FID); clf;
                set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
                for igen=1:ngens,
                    h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                        1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                    xlim([minx maxx]); ylim([miny maxy]); title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis on;
                    if strcmp(lrb,'r'), img = imgr{igen,ishift}; end
                    if strcmp(lrb,'l'), img = imgl{igen,ishift}; end
                    imagesc(xv,yv,img',[0 max_int]); drawnow; axis square;
                    %             patch([3*maxx/8 5*maxx/8 maxx/2 3*maxx/8], [3*maxx/8 3*maxx/8 5*maxx/8 3*maxx/8],'w');
                    plot([0 0], [0 1*maxx/8],'w');
                    plot([-1*maxx/8 1*maxx/8], [0 0],'w');
                    set(h,'FontSize',params.axisfontsize); set(gca,'FontSize',params.axisfontsize);
                    set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
                    set(gca,'YTick',xaxis); set(gca,'YTickLabel',xtick_lab);
                    %             if ~mod(4-1,params.k(1)), set(gca,'Visible','On'); else set(gca,'Visible','Off'); end
                    if ~mod(igen-1,params.k(2)), ylabel('distance [mm]','FontSize',labelfntsize); end
                    if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]','FontSize',labelfntsize); end
                    colormap('default'); colormap([0 0 0 ; jet]);
                end
                pos = get(gca,'Position');
                ch = colorbar('location','EastOutside'); set(get(ch,'YLabel'),'String','occurence [-]');
                set(ch,'Position',[0.93 pos(2) 0.01 pos(4)]);
                if params.courtship,
                    titl = 'Female Pos. ';
                else
                    titl = 'Fly Pos. ';
                end
                if shift<0,
                    titl = [titl 'BEFORE'];
                elseif shift>0,
                    titl = [titl 'AFTER'];
                else
                    titl = [titl 'WHILE'];
                end
                if params.courtship,
                    titl = [titl ' Male Ext. '];
                else
                    titl = [titl ' other Fly Ext. '];
                end
                if strcmp(lrb,'r'), titl = [titl 'RIGHT Wing']; else titl = [titl 'LEFT Wing']; end
                titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
                titl = [titl ', \omega > ' num2str(min_vphi/dts{1}(1),'%3.1f') '^\circ/s'];
                titl = [titl ', v_{orth} > ' num2str(min_vortho,'%3.1f') 'mm/s'];
                titl = [titl ', v_{para} > ' num2str(min_vparal,'%3.1f') 'mm/s'];
                set(FID,'Name',titl); mtit(titl,'FontSize',12,'yoff',.04);

                if params.pdf,
                    %             if (ilrb == 1) && (ishift == 1),
                    %                 print('-f130','-dpsc2',params.PSFileN);
                    %             else
                    print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN);
                    %             end
                end
            end
        end
    end
   
    if numel(find(params.plots.heatmaps.wingrel_vec == 2)),
        % PLOT HEATMAPS OF RELATIVE FLY POSTION WHILE EXTENDING A WING
        % DIFFERENCE BETWEEN LEFT AND RIGHT WING
        for ishift=2:2,
            shift = shift_vec(ishift);
            if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
            figure(FID); clf;
            col = jet(64); for j=1:2, col(ceil(length(col)/2)+j-1,:) = [0 0 0]; end
            set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
            for igen=1:ngens,
                h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                    1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
                xlim([minx maxx]); ylim([miny maxy]); title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis on;
                img = imgl{igen,ishift} - imgr{igen,ishift};
                imagesc(xv,yv,img',[-ceil(max_int/2) ceil(max_int/2)]); drawnow; axis square;
                plot([0 0], [0 1*maxx/8],'w');
                plot([-1*maxx/8 1*maxx/8], [0 0],'w');
                set(h,'FontSize',params.axisfontsize); set(gca,'FontSize',params.axisfontsize);
                set(gca,'XTick',xaxis); set(gca,'XTickLabel',xtick_lab);
                set(gca,'YTick',xaxis); set(gca,'YTickLabel',xtick_lab);
                if ~mod(igen-1,params.k(2)), ylabel('distance [mm]','FontSize',labelfntsize); end
                if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]','FontSize',labelfntsize); end
                colormap(col);
            end
            pos = get(gca,'Position');
            ch = colorbar('location','EastOutside'); set(get(ch,'YLabel'),'String','occurence [-]');
            set(ch,'Position',[0.93 pos(2) 0.01 pos(4)]);
            if params.courtship, titl = 'Female Pos. '; else titl = 'Fly Pos. '; end
            if shift<0,
                titl = [titl 'BEFORE'];
            elseif shift>0,
                titl = [titl 'AFTER'];
            else
                titl = [titl 'WHILE'];
            end
            if params.courtship, titl = [titl ' Male Ext. ']; else titl = [titl ' other Fly Ext. ']; end
            titl = [titl 'LEFT-RIGHT Wing'];
            titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
            titl = [titl ', \omega > ' num2str(min_vphi/dts{1}(1),'%3.1f') '^\circ/s'];
            titl = [titl ', v_{orth} > ' num2str(min_vortho,'%3.1f') 'mm/s'];
            titl = [titl ', v_{para} > ' num2str(min_vparal,'%3.1f') 'mm/s'];
            set(FID,'Name',titl); mtit(titl,'FontSize',12,'yoff',.04);
            if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
        end
    end

    if numel(find(params.plots.heatmaps.wingrel_vec == 3)),
        % PLOT CIRCLING DIRECTION DURING WING EXTENSION PHASES
        dw = 0.1667; sh = 2*dw;
        data = rotdir1 + rotdir2;
        ishift = 2;
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for ilrb=1:2,
            dp = sign(ilrb-1.5)*dw;
            for igen=1:ngens,
                patch([igen+dp-sh,igen-sh,igen-sh,igen+dp-sh,igen+dp-sh],...
                    [0,0,data(ishift,ilrb,igen,1),data(ishift,ilrb,igen,1),0],params.flycol.win);
                patch([igen+dp-sh,igen-sh,igen-sh,igen+dp-sh,igen+dp-sh],...
                    [0,0,data(ishift,ilrb,igen,2),data(ishift,ilrb,igen,2),0],[.71 .82 .49]);
                hold on;
            end
        end
        text(1-dw-sh+0.03,0.1,'L'); text(1-sh+0.03,0.1,'R');
        h = legend('cclkw','clkw'); legend(h,'boxoff'); axis([dw ngens+dw -1 1]);
        set(gca,'FontSize',params.axisfontsize);
        set(gca,'XTick',dw:1:ngens+dw); set_xtick_label(gen_ident,90,params.cross_to,params.axisfontsize);
        set(gca,'YTick',-1:.5:1); set(gca,'YTickLabel',{100,50,0,50,100});
        set(gca,'YGrid','on'); ylabel('clockwise <-- fly pairs [%] --> counter-clockwise')
        titl = 'Circling Direction While Wing Extension';
        titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
        titl = [titl ', \omega > ' num2str(min_vphi/dts{1}(1),'%3.1f') '^\circ/s'];
        titl = [titl ', v_{orth} > ' num2str(min_vortho,'%3.1f') 'mm/s'];
        titl = [titl ', v_{para} > ' num2str(min_vparal,'%3.1f') 'mm/s'];
        set(FID,'Name',titl); mtit(titl,'FontSize',12,'yoff',.04);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

% HEATMAP TIME SERIES OF RELATIVE FLY ORIENTATION DURING WING-EXTENSION
bool_test = 0;
if bool_test,
    extfli = 'ext';
    maxx = 10; % chamber.width, Range of distances of interest
    max_int = 200; % max. frequency in histogram
    add = 30; % add. n frames before/after each sequence
    % phi0.min = 40; phi0.max = 60; % borders for wing angle (rel. to major body-axis)
    % phi0.min = 60; phi0.max = 100; % borders for wing angle (rel. to major body-axis)
    phi0.min = 70; phi0.max = 90; % borders for wing angle (rel. to major body-axis)
    tstep = 5; % divide an event into tstep time steps
    lrb_vec = ['l' 'r'];

    for ilrb=1:2,
        lrb = lrb_vec(ilrb);
        fn = [path addname num2str(ngens) '_gens_wing_' extfli '_' lrb '_add_' ...
            num2str(add) '_minmaxphi_' num2str(phi0.min) '_' num2str(phi0.max) '_frames'];
        fn = [fn '_ob_orient'];
        ngens = max(genotypes);
        dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
        o1 = cell(1,ngens); o2 = o1;
        fid = fopen([fn '.mat'],'r');
        if (fid < 0) || params.analyze_new,
            h1 = waitbar(0,'Analyzing Data');
            for igen=1:ngens,
                waitbar(igen/ngens,h1,gen_ident{igen});
                [ob1,ob2] = extract_wing_orient(igen,wings,extfli,lrb,round(add/dt),0,phi0,genotypes,sorted_names,copu);
                o1{igen} = ob1; o2{igen} = ob2;
            end
            close(h1);
            save(fn,'o1','o2');
        else
            load(fn);
        end

        % RELATIVE FLY POSITION HEATMAPS FOR WING EXTENSION PHASES
        % OVER TIME
        nt = add/tstep*2;
        minx = -maxx; miny = -maxx; maxy = maxx; dv = maxx/15;
        xv = minx:dv:maxx; yv = miny:dv:maxy;
        ngens = max(genotypes); %k = factor(ngens);
        img = cell(ngens,nt); phi_velo1 = cell(ngens,2); phi_velo2 = phi_velo1;
        for igen=1:ngens,
            if ~params.courtship || params.oneobj,
                [img,phi_velo1] = analyze_wing_mov(o1,phi_velo1,img,igen,add,tstep,min_vortho,min_vparal,min_vphi,xv,yv,chamber,params.oneobj);
            end
            if ~params.courtship || ~params.oneobj,
                [img,phi_velo2] = analyze_wing_mov(o2,phi_velo2,img,igen,add,tstep,min_vortho,min_vparal,min_vphi,xv,yv,chamber,params.oneobj);
            end
            if strcmp(lrb,'r'), imgr = img; end
            if strcmp(lrb,'l'), imgl = img; end
        end
    end
    
    % PLOT HEATMAP TIME SERIES
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    col = jet; %for j=1:2, col(ceil(length(col)/2)+j-1,:) = [0 0 0]; end
    for l1=1:nt,
        for igen=1:ngens,
            dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            xlim([minx maxx]); ylim([miny maxy]); title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize); axis on;
            img = imgl{igen,l1} - imgr{igen,l1};
            imagesc(xv,yv,img',[-ceil(max_int/4) ceil(max_int/4)]); drawnow; axis square;
            plot([0 0], [0 1*maxx/8],'w');
            plot([-1*maxx/8 1*maxx/8], [0 0],'w');
            set(h,'FontSize',params.axisfontsize);
            colormap(col);
            title(['t = ' num2str((-add+((l1-1)*tstep+add/2))*dt) ' s']);
        end
        if ~mod(igen-1,params.k(2)), ylabel('distance [mm]'); end
        if (igen > (params.k(1)-1)*params.k(2)), xlabel('distance [mm]'); end
        pos = get(gca,'Position');
        ch = colorbar('location','EastOutside'); set(get(ch,'YLabel'),'String','occurence [-]');
        set(ch,'Position',[0.93 pos(2) 0.01 pos(4)]);
        titl = 'Histograms of Fly Position';
        if strcmp(lrb,'r'), titl = [titl ' RIGHT Wing Extended']; else titl = [titl ' LEFT Wing Extended']; end
        titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
        titl = [titl ', \omega > ' num2str(min_vphi/dts{1}(1),'%3.1f') '^\circ/s'];
        titl = [titl ', v_{orth} > ' num2str(min_vortho,'%3.1f') 'mm/s'];
        titl = [titl ', v_{para} > ' num2str(min_vparal,'%3.1f') 'mm/s'];
        set(FID,'Name',titl); mtit(titl,'FontSize',14','yoff',.04);
        if params.pdf,
            %             if (ilrb == 1) && (l1 == 1),
            %                 print('-f130','-dpsc2',params.PSFileN);
            %             else
            print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN);
            %             end
        end
    end
end


% FLY DISTANCE OVER TIME FOR WING EXTENSION PHASES
if ~params.oneobj && params.plots.stat.wingextflydist,
    extfli = 'ext';
    add = 1; %[s]; add. n frames before/after each sequence
    % phi0.min = 40; phi0.max = 60; % borders for wing angle (rel. to major body-axis)
    % phi0.min = 60; phi0.max = 100; % borders for wing angle (rel. to major body-axis)
    phi0.min = 70; phi0.max = 90; % borders for wing angle (rel. to major body-axis)
    lrb_vec = ['l' 'r'];
    shift = 0;
    ngens = max(genotypes);

    pred_mean = zeros(2,ngens); postd_mean = pred_mean;
    dis_vec = cell(2,ngens); dur_vec = dis_vec;
    for ilrb=1:2,
        lrb = lrb_vec(ilrb);
        % EXTRACT RELATIV ORIENTATION AND POSITION OF ONE FLY TOWARDS
        % THE OTHER FLY, WHILE THE OTHER FLY IS EXTENDING A WING
        fn = [path addname num2str(ngens) '_gens_wing_' extfli '_' lrb '_add_' ...
            num2str(add) '_minmaxphi_' num2str(phi0.min) '_' num2str(phi0.max) '_frames'];
        if shift<0,
            fn = [fn '_before_ob_orient'];
        elseif shift>0
            fn = [fn '_after_ob_orient'];
        else
            fn = [fn '_ob_orient'];
        end
        o1 = cell(1,ngens); o2 = o1;
        fid = fopen([fn '.mat'],'r');
        if (fid < 0) || params.analyze_new,
            h1 = waitbar(0,'Analyzing Data');
            for igen=1:ngens,
                waitbar(igen/ngens,h1,gen_ident{igen});
                dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
                [ob1,ob2] = extract_wing_orient(igen,wings,extfli,lrb,round(add/dt),shift,phi0,genotypes,sorted_names,copu);
                o1{igen} = ob1; o2{igen} = ob2;
            end
            close(h1);
            save(fn,'o1','o2');
        else
            load(fn);
        end

        mint = -add; maxt = add;
        col = jet(15);
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for igen=1:ngens,
            dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
            xt = -add:dt:add; if mod(length(xt),2), xt = [xt add+dt]; end;
            add1 = round(length(xt)/2);
            h = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            cnt = 0;
            dis_arr = []; dis_vec{ilrb,igen} = []; dur_vec{ilrb,igen} = [];
            mdis = zeros(1,length(xt)); sdis = mdis; rmdis = mdis; rsdis = mdis;
            % Fly 1
            if ~params.courtship || params.oneobj,
                for j=1:length(o1{igen}.tim),
                    if numel(o1{igen}.tim{j}),
                        for l=1:length(o1{igen}.tim{j}),
                            if length(o1{igen}.dis{j}{l})>2*add1,
                                cnt = cnt + 1;
                                % In case of a single fly analyze with respect to arena
                                % center
                                if params.oneobj,
                                    r = sqrt((o1{igen}.x1{j}{l}(1:length(xt))-chamber.width/2).^2 + (o1{igen}.y1{j}{l}(1:length(xt))-chamber.height/2).^2);
                                else
                                    r = o1{igen}.dis{j}{l}(1:length(xt));
                                end
                                dis = conv2(r,[.25 .5 .25],'same');
                                % Normalize to center frame
                                dis_arr = [dis_arr ; dis / dis(add1)];
                                % Fly distance and time at center frame
                                dis_vec{ilrb,igen} = [dis_vec{ilrb,igen} o1{igen}.dis{j}{l}(add1+1)];
                                dur_vec{ilrb,igen} = [dur_vec{ilrb,igen} o1{igen}.tim{j}{l}(end-add1)-o1{igen}.tim{j}{l}(add1+1)];
                            end
                        end
                    end
                end
            end
            % Fly 2
            if ~params.courtship || ~params.oneobj,
                for j=1:length(o2{igen}.tim),
                    if numel(o2{igen}.tim{j}),
                        for l=1:length(o2{igen}.tim{j}),
                            if length(o2{igen}.dis{j}{l})>2*add1,
                                cnt = cnt + 1;
                                % In case of a single fly analyze with respect to arena
                                % center
                                if params.oneobj,
                                    r = sqrt(o2{igen}.x1{j}{l}(1:length(xt)).^2 + o2{igen}.y1{j}{l}(1:length(xt)).^2);
                                else
                                    r = o2{igen}.dis{j}{l}(1:length(xt));
                                end
                                dis = conv2(r,[.25 .5 .25],'same');
                                % Normalize to center frame
                                dis_arr = [dis_arr ; dis / dis(add1)];
                                % Fly distance and time at center frame
                                dis_vec{ilrb,igen} = [dis_vec{ilrb,igen} o2{igen}.dis{j}{l}(add1+1)];
                                dur_vec{ilrb,igen} = [dur_vec{ilrb,igen} o2{igen}.tim{j}{l}(end-add1)-o2{igen}.tim{j}{l}(add1+1)];
                            end
                        end
                    end
                end
            end
            % Statistics
            mdis = mean(dis_arr,1); sdis = std(dis_arr,[],1) / sqrt(cnt);
            q25dis = quantile(dis_arr,.25,1); q75dis = quantile(dis_arr,.75,1);
            % Statistics before/after wing extension
            if length(mdis) > add1,
                pred_mean(ilrb,igen) = mean(mdis(1:add1-1));
                postd_mean(ilrb,igen) = mean(mdis(add1+1:end));
            end
            
            % PLOT FLY DISTANCE OVER TIME
            randvec = randperm(length(dists{igen}.distc)-2*add1);
            lab = [wings.ext.l{igen}.obj1.tim wings.ext.r{igen}.obj1.tim wings.ext.b{igen}.obj1.tim];
            nrand = 100; rdis_arr = []; cnt1 = 0; l = 0;
            while (cnt1 < nrand) && (l < length(randvec)),
                l = l + 1;
                cnt1 = cnt1 + 1;
                st = randvec(l);
                rdis_arr = [rdis_arr ; dists{igen}.distc(st:st+2*add1-1) / dists{igen}.distc(st+add1-1)];
            end
            rmdis = mean(rdis_arr,1); rsdis = std(rdis_arr,[],1) / sqrt(nrand);
            q25rdis = quantile(rdis_arr,.25,1); q75rdis = quantile(rdis_arr,.75,1);
            if length(xt) == length(mdis),
                fill([xt -xt],[mdis mdis(end:-1:1)+sdis(end:-1:1)],[0.4 1 0.4],'EdgeColor','none'); hold on;
                fill([xt -xt],[mdis mdis(end:-1:1)-sdis(end:-1:1)],[0.4 1 0.4],'EdgeColor','none'); hold on;
                plot(xt,mdis,'LineWidth',4,'Color',[0 .5 0]); hold on;
            end
            set(h,'FontSize',params.axisfontsize);
            ylim([.5 2]);
            if ~mod(igen-1,params.k(2)), ylabel('average normalized distance'); end
            if igen>2*params.k(2), xlabel('time [s]'); end
            %         if igen==ngens, legend('ext','no ext'); end
            axis([mint maxt .5 2]);
            text(mint+(maxt-mint)*0.1,0.1,['n=' num2str(cnt)]);
        end
        titl = 'Distance of Flies under ';
        if strcmp(lrb,'r'), titl = [titl 'RIGHT Wing Extension']; else titl = [titl 'LEFT Wing Extension']; end
        titl = [titl ' (mean, sem)'];
        titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
        set(FID,'Name',titl); mtit(titl,'FontSize',14,'yoff',.0);
        if params.pdf, print(['-f' num2str(FID)] ,'-dpsc2',appnd,params.PSFileN); end
    end

    % PLOT NORMALIZED DISTANCE BEFORE/AFTER WING EXTENSION
    pred_mean = mean(pred_mean,1); postd_mean = mean(postd_mean,1);
    normd_ratio = postd_mean - pred_mean;
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    bar([pred_mean' postd_mean']-1); colormap([params.flycol.win ; params.flycol.los]); hold on;
    axis([0.5 ngens+.5 -1 1]);
    set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,90,params.cross_to,params.axisfontsize);
    set(gca,'FontSize',params.axisfontsize);
    ylabel('normalized distance [-]'); h = legend('prior','post'); legend(h,'boxoff');
    titl = 'Mean Norm. Distance - Prior & Post Wing Extension';
    mtit(titl,'FontSize',params.axisfontsize,'yoff',.0); set(FID,'Name',titl);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end

    % PLOT DIFFERENCE OF NORMALIZED DISTANCE 
    % (AFTER - BEFORE) WING EXTENSION
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    bar(normd_ratio'); colormap([params.flycol.win ; params.flycol.los]); hold on;
    axis([0.5 ngens+.5 -1 1]); set(gca,'XTick',1:ngens); set_xtick_label(gen_ident,90,params.cross_to,params.axisfontsize);
    set(gca,'FontSize',params.axisfontsize);
    ylabel('normalized distance [-]');
    titl = 'Difference of Norm. Distance - Post/Prior Wing Extension';
    mtit(titl,'FontSize',params.axisfontsize,'yoff',.0); set(FID,'Name',titl);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end

    % % PLOT DISTANCE AT WING EXTENSION OVER TIME WHEN IT OCCURRED 
    % for ilrb=1:2,
    %     lrb = lrb_vec(ilrb);
    %     figure(FID); clf;
    %     set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    %     for igen=1:ngens,
    %         subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
    %                             1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
    %         title(cell2mat(gen_ident(igen)));
    %         plot(dis_vec{ilrb,igen},dur_vec{ilrb,igen},'.','MarkerSize',6,'Color',[0 0 .8]);
    %         axis([0 50 0 5]); grid on;
    %         text(25,4.5,['n=' num2str(numel(dur_vec{ilrb,igen}))]);
    %         if ~mod(igen-1,params.k(1)), xlabel('fly distance [mm]'); end
    %         if ~mod(igen-1,params.k(2)), ylabel('duration [s]'); end
    %     end
    %     titl = 'Fly Distance at ';
    %     if strcmp(lrb,'r'), titl = [titl 'RIGHT ']; else titl = [titl 'LEFT ']; end
    %     titl = [titl 'Wing Ext. <-> Wing Ext. Duration'];
    %     set(FID,'Name',titl); mtit(titl,'FontSize',14','yoff',.04);
    %     if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    % end
end

% FLY ORIENTATION POLAR PLOTS FOR WING EXTENSION PHASES
bool_test = 0;
if bool_test,
    extfli = 'ext';
    add = 30; %[s]; add. n frames before/after each sequence
    sh = 5;
    lrb = 'l';
    shift = 0;
    % phi0.min = 40; phi0.max = 60; % borders for wing angle (rel. to major body-axis)
    % phi0.min = 60; phi0.max = 100; % borders for wing angle (rel. to major body-axis)
    ngens = max(genotypes);
    maxr = 5; % chamber.width, Range of distances of interest

    for ilrb=1:2,
        lrb = lrb_vec(ilrb);
        fn = [path addname num2str(ngens) '_gens_wing_' extfli '_' lrb '_add_' ...
            num2str(add) '_minmaxphi_' num2str(phi0.min) '_' num2str(phi0.max) '_frames'];
        if shift<0,
            fn = [fn '_before_ob_orient'];
        elseif shift>0
            fn = [fn '_after_ob_orient'];
        else
            fn = [fn '_ob_orient'];
        end
        o1 = cell(1,ngens); o2 = o1;
        fid = fopen([fn '.mat'],'r');
        if (fid < 0) || params.analyze_new,
            h1 = waitbar(0,'Analyzing Data');
            for igen=1:ngens,
                waitbar(igen/ngens,h1,gen_ident{igen});
                dt = sum(dts{igen}.*nframes{igen})/sum(nframes{igen});
                [ob1,ob2] = extract_wing_orient(igen,wings,extfli,lrb,round(add/dt),shift,phi0,genotypes,sorted_names,copu);
                o1{igen} = ob1; o2{igen} = ob2;
            end
            close(h1);
            save(fn,'o1','o2');
        else
            load(fn);
        end

        col = jet(15);
        if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
        figure(FID); clf;
        set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
        for igen=1:ngens,
            h1 = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
                1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
            dori = []; r = [];
            if ~params.courtship || params.oneobj,
                for j=1:length(o1{igen}.tim),
                    if numel(o1{igen}.tim{j}),
                        for l=1:length(o1{igen}.tim{j}),
                            if length(o1{igen}.dis{j}{l})>2*add,
                                if params.oneobj,
                                    phi = atan2(o1{igen}.y1{j}{l}(add+sh),o1{igen}.x1{j}{l}(add+sh)) * 180/pi;
                                    r = [r ; sqrt(o1{igen}.x1{j}{l}(add+sh).^2 + o1{igen}.y1{j}{l}(add+sh).^2)];
                                else
                                    phi = atan2((o1{igen}.y2{j}{l}(add+sh)-o1{igen}.y1{j}{l}(add+sh)),(o1{igen}.x2{j}{l}(add+sh)-o1{igen}.x1{j}{l}(add+sh))) * 180/pi;
                                    r = [r ; o1{igen}.dis{j}{l}(add+sh)];
                                end
                                dori = [dori ; phi - o1{igen}.do1{j}{l}(add+sh) + 90];
                            end
                        end
                    end
                end
            end
            if ~params.courtship || ~params.oneobj,
                for j=1:length(o2{igen}.tim),
                    if numel(o2{igen}.tim{j}),
                        for l=1:length(o2{igen}.tim{j}),
                            if length(o2{igen}.dis{j}{l})>2*add,
                                phi = atan2((o2{igen}.y2{j}{l}(add+sh)-o2{igen}.y1{j}{l}(add+sh)),(o2{igen}.x2{j}{l}(add+sh)-o2{igen}.x1{j}{l}(add+sh))) * 180/pi;
                                r = [r ; o2{igen}.dis{j}{l}(add+sh)];
                                dori = [dori ; phi - o2{igen}.do1{j}{l}(add+sh) + 90];
                                %                             dori = [dori ; o2{igen}.do2{j}{l}(add+sh) - o2{igen}.do1{j}{l}(add+sh)];
                                %                             r = [r ; o2{igen}.dis{j}{l}(add+sh)];
                            end
                        end
                    end
                end
            end
            in = find(dori > 180); if numel(in), dori(in) = dori(in) - 360; end
            in = find(dori < -180); if numel(in),dori(in) = dori(in) + 360; end
            dori = dori * pi/180;

            in = find(r <= maxr); r = r(in); dori = dori(in);
            mean_ori = atan2(mean(sin(dori)),mean(cos(dori)));
            stdv_ori = atan2(std(sin(dori)),std(cos(dori)));
            h = polar(dori,r,'.'); hold on; axis square;
            set(h,'Color',params.flycol.win,'LineWidth',2);
            h = polar([mean_ori mean_ori],[0 maxr]);
            set(h,'Color',[0 0 0.8],'LineWidth',4);
            h = polar([mean_ori-stdv_ori mean_ori-stdv_ori],[0 maxr],'.');
            set(h,'Color','r','LineWidth',2,'LineStyle','-');
            h = polar([mean_ori+stdv_ori mean_ori+stdv_ori],[0 maxr],'.');
            set(h,'Color','r','LineWidth',2,'LineStyle','-');
            set(h1,'FontSize',params.axisfontsize);
        end
        titl = 'Mean Orientation of Flies ';
        if (i == 1), titl = [titl 'Before ']; elseif (i == 2), titl = [titl 'During ']; else titl = [titl 'After ']; end
        if strcmp(lrb,'r'), titl = [titl 'Right ']; else titl = [titl 'Left ']; end
        titl = [titl 'Wing Extension'];
        titl = [titl ' >' num2str(phi0.min) '^o <' num2str(phi0.max) '^o'];
        set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize,'yoff',.0);
        if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
    end
end

%% HISTOGRAMS OF WING EXTENSION DURATIONS
if (params.plots.stat.leftwing || params.plots.stat.rightwing || ...
    params.plots.stat.onewing) && ~params.plots.movieclips && ...
    numel(find(params.plots.stat.wing_vec == 2)),
    mint = 0; maxt = params.max_frames/60; dt1 = 2;
    t = mint:dt1:maxt; xtick_lab = cell(1,length(t));
    for i=1:length(t), xtick_lab{i} = num2str(t(i)); end
    th = mint:dt1:maxt;
    if ~FID, FID = 1; appnd = []; else FID = FID + 1; appnd = '-append'; end
    figure(FID); clf;
    set(FID,'PaperOrientation','landscape','PaperPositionMode','manual','PaperPosition',[0 0 11 8.5]);
    h = cell(1,ngens); maxy = 0;
    for igen=1:ngens,
        nmov = length(wings.ext.l{igen}.obj1.number);
        h2d0 = zeros(nmov,length(th));
        ci_1.l = [0 cumsum(wings.ext.l{igen}.obj1.number)];
        ci_1.r = [0 cumsum(wings.ext.r{igen}.obj1.number)];
        ci_2.l = [0 cumsum(wings.ext.l{igen}.obj2.number)];
        ci_2.r = [0 cumsum(wings.ext.r{igen}.obj2.number)];
        for j=1:nmov,
            % COLLECT WING EXTENSION PERIODS FOR EACH MOVIE
            dt_1.l = ones(1,wings.ext.l{igen}.obj1.number(j))*dts{igen}(j);
            dt_2.l = ones(1,wings.ext.l{igen}.obj2.number(j))*dts{igen}(j);
            dt_1.r = ones(1,wings.ext.r{igen}.obj1.number(j))*dts{igen}(j);
            dt_2.r = ones(1,wings.ext.r{igen}.obj2.number(j))*dts{igen}(j);
            if ~(numel(dt_2.l)+numel(dt_2.r)),
                data = [wings.ext.l{igen}.obj1.len(ci_1.l(j)+1:ci_1.l(j+1)).*dt_1.l ...
                    wings.ext.r{igen}.obj1.len(ci_1.r(j)+1:ci_1.r(j+1)).*dt_1.r];
            else
                data = [wings.ext.l{igen}.obj1.len(ci_1.l(j)+1:ci_1.l(j+1)).*dt_1.l ...
                    wings.ext.l{igen}.obj2.len(ci_2.l(j)+1:ci_2.l(j+1)).*dt_2.l ...
                    wings.ext.r{igen}.obj1.len(ci_1.r(j)+1:ci_1.r(j+1)).*dt_1.r ...
                    wings.ext.r{igen}.obj2.len(ci_2.r(j)+1:ci_2.r(j+1)).*dt_2.r];
            end
            % Histogram of wing extension periods
            h2d0(j,:) = hist(data,th);
        end
        % MEAN HISTOGRAM FOR EACH GENOTYPE
        h{igen} = mean(h2d0,1);
        if max(h{igen}) > maxy, maxy = ceil(max(h{igen})); end
    end
    % PLOT HISTOGRAMS
    if ~params.autoscale, maxy = 100; end
    for igen=1:ngens,
        h1 = subplot('Position',[mod(igen-1,params.k(2))/(params.k(2)*1.1)+.05, 0.9-fix((igen-1)/params.k(2))/(params.k(1)*1.1)-1/(params.k(1)*1.5), ...
            1/(params.k(2)*1.5), 1/(params.k(1)*1.5)]); hold on;
        title(cell2mat(gen_ident(igen)),'FontSize',params.axisfontsize);
        bar(th+dt1/2,h{igen},params.barwidth); colormap(params.flycol.win);
        %     colormap([0 0 0 ; 0 0 1]);
        set(h1,'FontSize',params.axisfontsize);
        axis([mint maxt 0 maxy]); set(gca,'XTick',mint:dt1:maxt);
        set(gca,'XTick',mint:dt1:maxt); set(gca,'XTickLabel',xtick_lab);
        xlim([mint maxt]);
        if (igen > (params.k(1)-1)*params.k(2)), xlabel('time [s]'); end
        if ~mod(igen-1,params.k(2)), ylabel('frequency of wing extensions'); end
    end
    titl = 'MEAN Wing Extension Frequency';
    set(FID,'Name',titl); mtit(titl,'FontSize',params.axisfontsize+2,'yoff',.04);
    if params.pdf, print(['-f' num2str(FID)],'-dpsc2',appnd,params.PSFileN); end
end

end