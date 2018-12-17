function [mov,tcnt] = write_sequence_crop(obj,obj1,obj2,roi,tcnt,params,scale,min_dist,framen,fig,img,dt,mov)

fac = 0.5; % scaling factor for resolution
center_over = 1; % frame center over 1 = fly 1, 2 = fly 2, 0 = mean

flywidth = .5;

img = imresize(img,fac);
scale.x = scale.x / fac;
scale.y = scale.y / fac;
roi.cols(1) = round(roi.cols(1) * fac);
roi.rows(1) = round(roi.rows(1) * fac);
r = 128/2; %ceil(min_dist/scale.x);
%% Crop area around center of both flies
if obj.i1112.xc < 1, obj.i1112.xc = obj.i1111.xc; end
if obj.i1112.yc < 1, obj.i1112.yc = obj.i1111.yc; end
if center_over == 1,
    xm = obj.i1111.xc;
    ym = obj.i1111.yc;    
elseif center_over == 2,
    xm = obj.i1112.xc;
    ym = obj.i1112.yc;    
else
    xm = (obj.i1111.xc + obj.i1112.xc) / 2;
    ym = (obj.i1111.yc + obj.i1112.yc) / 2;
end

minc = floor(xm/scale.x) - r + roi.cols(1); 
maxc = floor(xm/scale.x) + r + roi.cols(1)-1; 
minr = floor(ym/scale.y) - r + roi.rows(1); 
maxr = floor(ym/scale.y) + r + roi.rows(1)-1; 

img0 = img(max(minr,1):min(maxr,size(img,1)),...
          max(minc,1):min(maxc,size(img,2)),:);

if minc < 1, 
    addcr.cmin = 1 - minc;
else
    addcr.cmin = 0; 
end
if minr < 1, 
    addcr.rmin = 1 - minr;
else
    addcr.rmin = 0;
end
if maxc > size(img,2), 
    addcr.cmax = maxc - size(img,2); 
else
    addcr.cmax = 0; 
end
if maxr > size(img,1), 
    addcr.rmax = maxr - size(img,1); 
else
    addcr.rmax = 0; 
end
img1 = ones(addcr.rmin+addcr.rmax+size(img0,1),...
             addcr.cmin+addcr.cmax+size(img0,2),3);
img1(addcr.rmin+1:2*r-addcr.rmax,...
             addcr.cmin+1:2*r-addcr.cmax,:) = img0;
       
figure(fig);
set(fig,'DoubleBuffer','on');
set(gca,'xlim',[0 2*r],'ylim',[0 2*r], 'NextPlot','replace','Visible','off')
img1(img1 < 0) = 0; img1(img1 > 1) = 1;
image(img1);
hold on; axis equal off;

sec  = floor((framen-2)*dt);
tho  = floor(((framen-2)*dt - sec)*1000);
text(2*r-25,5,[num2str(sec,'%04d') '.' num2str(tho,'%03d') ' s'],'Color','k');

minr = max(minr,1) - addcr.rmin - 1;
minc = max(minc,1) - addcr.cmin - 1;      

%%
% The first two _line_ commands draw small grey and
% yellow rectangles to mark the fly locations on the
% screen. The next two _line_ commands draw the blue
% and red triangles to mark the fly directions inside
% each chambers.

line( [ obj.i1111.xc - 0.1 , obj.i1111.xc - 0.1 , ...
    obj.i1111.xc + 0.1 , obj.i1111.xc + 0.1 , ...
    obj.i1111.xc - 0.1 ]/ scale.x + roi.cols(1) - minc, ...
    [ obj.i1111.yc - 0.1 , obj.i1111.yc + 0.1 , ...
    obj.i1111.yc + 0.1 , obj.i1111.yc - 0.1 , ...
    obj.i1111.yc - 0.1 ]/ scale.y + roi.rows(1) - minr, ...
    'Color',[0.2 0.2 0.2]);
line( [ obj.i1112.xc - 0.1 , obj.i1112.xc - 0.1 , ...
    obj.i1112.xc + 0.1 , obj.i1112.xc + 0.1 , ...
    obj.i1112.xc - 0.1 ]/ scale.x + roi.cols(1) - minc, ...
    [ obj.i1112.yc - 0.1 , obj.i1112.yc + 0.1 , ...
    obj.i1112.yc + 0.1 , obj.i1112.yc - 0.1 , ...
    obj.i1112.yc - 0.1] / scale.y + roi.rows(1) - minr, ...
    'Color', [.6 .6 0] );

fc1 = flywidth*(cos(obj.i1111.head+params.hpi));
fs1 = flywidth*(sin(obj.i1111.head+params.hpi));
fc2 = flywidth*(cos(obj.i1112.head+params.hpi));
fs2 = flywidth*(sin(obj.i1112.head+params.hpi));
line( [ obj.i1111.xh , obj.i1111.xt - fc1 , ...
    obj.i1111.xt + fc1 , obj.i1111.xh ] / ...
    scale.x + roi.cols(1) - minc, ...
    [ obj.i1111.yh , obj.i1111.yt - fs1 , ...
    obj.i1111.yt + fs1 , obj.i1111.yh ] / ...
    scale.y + roi.rows(1) - minr, ...
    'Color', [.2 .2 1], 'LineWidth', 2 );
line( [ obj.i1112.xh , obj.i1112.xt - fc2, ...
    obj.i1112.xt + fc2 , obj.i1112.xh ] / ...
    scale.x + roi.cols(1) - minc, ...
    [ obj.i1112.yh , obj.i1112.yt - fs2, ...
    obj.i1112.yt + fs2 , obj.i1112.yh] / ...
    scale.y + roi.rows(1) - minr, ...
    'Color', [1 .2 .2], 'LineWidth', 2 );

%% Wings - draw line from fly body center to wing tip as
%% measured
f1rx = cos(obj.i1111.head-obj.i1111.phir+pi)*obj.i1111.r+obj.i1111.xc;
f1ry = sin(obj.i1111.head-obj.i1111.phir+pi)*obj.i1111.r+obj.i1111.yc;
f2rx = cos(obj.i1112.head-obj.i1112.phir+pi)*obj.i1112.r+obj.i1112.xc;
f2ry = sin(obj.i1112.head-obj.i1112.phir+pi)*obj.i1112.r+obj.i1112.yc;
f1lx = cos(obj.i1111.head+obj.i1111.phil+pi)*obj.i1111.l+obj.i1111.xc;
f1ly = sin(obj.i1111.head+obj.i1111.phil+pi)*obj.i1111.l+obj.i1111.yc;
f2lx = cos(obj.i1112.head+obj.i1112.phil+pi)*obj.i1112.l+obj.i1112.xc;
f2ly = sin(obj.i1112.head+obj.i1112.phil+pi)*obj.i1112.l+obj.i1112.yc;
line( [ obj.i1111.xc , f1rx] / scale.x + roi.cols(1) - minc, ...
    [ obj.i1111.yc , f1ry] / scale.y + roi.rows(1) - minr, ...
    'Color', [.2 .9 .2], 'LineWidth', 2 );
line( [ obj.i1111.xc , f1lx] / scale.x + roi.cols(1) - minc, ...
    [ obj.i1111.yc , f1ly] / scale.y + roi.rows(1) - minr, ...
    'Color', [.9 .9 .2], 'LineWidth', 2 );
line( [ obj.i1112.xc , f2rx] / scale.x + roi.cols(1) - minc, ...
    [ obj.i1112.yc , f2ry] / scale.y + roi.rows(1) - minr, ...
    'Color', [.2 .9 .2], 'LineWidth', 2 );
line( [ obj.i1112.xc , f2lx] / scale.x + roi.cols(1) - minc, ...
    [ obj.i1112.yc , f2ly] / scale.y + roi.rows(1) - minr, ...
    'Color', [.9 .9 .2], 'LineWidth', 2 );
%%
% The following motion lines can only be drawn after
% more than seven frames (measurements) are made, for
% the obvious reason discussed above: the first
% acceleration vector is availabe only after six frames
% are processed.

if (framen > 7),
    
    for b=1:tcnt,
        if (tcnt > 10),
            c='r';
        else
            c='g';
        end
        if (b > 10),
            line( [ obj1.pos_x(framen-3+b-tcnt), ...
                obj1.pos_x(framen-2+b-tcnt) ] / ...
                scale.x + roi.cols(1) - minc, ...
                [ obj1.pos_y(framen-3+b-tcnt), ...
                obj1.pos_y(framen-2+b-tcnt) ] / ...
                scale.y + roi.rows(1) - minr, 'Color', 'g' );
            line( [ obj2.pos_x(framen-3+b-tcnt), ...
                obj2.pos_x(framen-2+b-tcnt) ] / ...
                scale.x + roi.cols(1) - minc, ...
                [ obj2.pos_y(framen-3+b-tcnt), ...
                obj2.pos_y(framen-2+b-tcnt) ] / ...
                scale.y + roi.rows(1) - minr, 'Color', 'g' );
        else
            line( [ obj1.pos_x(framen-3+b-tcnt), ...
                obj1.pos_x(framen-2+b-tcnt)] / ...
                scale.x + roi.cols(1) - minc, ...
                [ obj1.pos_y(framen-3+b-tcnt), ...
                obj1.pos_y(framen-2+b-tcnt)] / ...
                scale.y + roi.rows(1) - minr, 'Color', c );
            line( [ obj2.pos_x(framen-3+b-tcnt), ...
                obj2.pos_x(framen-2+b-tcnt)] / ...
                scale.x + roi.cols(1) - minc, ...
                [ obj2.pos_y(framen-3+b-tcnt), ...
                obj2.pos_y(framen-2+b-tcnt)] / ...
                scale.y + roi.rows(1) - minr, 'Color', c );
        end
    end
    tcnt = tcnt + 1;
    if (tcnt > 30),
        tcnt = 28;
    end
    
end
drawnow;

F = getframe(gca);
mov = addframe(mov,F);

end
