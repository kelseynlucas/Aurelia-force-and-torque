%July 12, 2018

%Reads the pressure fields and calculates force per unit depth, and
%torque.





%Will need to add interparc.m and polygeom.m to the Matlab file path before 
%running this script.

%Optionally need cbrewer.m if using its color palette (commented out by default)




clear all
close all


% %If you want a nicer palette for pressure fields, load colorBrewer red to 
% %blue diverging palette.  Make sure cbrewer.m is in the matlab path.
% mycolors=cbrewer('div','RdBu',128);                       
% %Flips palette so that red = high values
% fci = 0;
% for row = 0:1:127
%     fci = fci+1;
%     flippedcolors(fci,:) = mycolors(128-row,:);
% end






%%%%%START MANUAL PARAMETER SETTINGS%%%%%



%Set metadata: filepath information, scale factors, numberformat in the filenames,
%delimiters, and increment between frames used for pressure code.

numBoundPts = 100;      %Set the number of points you want to calculate force at


vecstub = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection/B';    %File path and name stub for vector fields
imstub = 'G:/Jellyfish turns/Aurelia turn/July 2018/Nov 17_2011_Aurelia 6cm_03/piv inc10/bitmap grey1/B';     %File path and name stub for images
blankstub = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection/newblank';      %File path and name stub for blanking boundaries
pressstub = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection/pressure_knl_';    %File path and name stub for smoothed pressure field file
pixscale = 0.124809;   %scale factor in mm/pixel
increment = 10;  %increment between file numbers
first = 1;  %number of first frame
last = 261;  %number of last frame
savePath = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection';       %set folder where data should be saved
imfileformat = 'bmp';

%After doing the rest of the settings,
%Set all three of these to "yes".
%Turn on the "return" in the figure(2) section. Run script. Answer the
%following based on the figures that appear.
%Is the figure(1) jellyfish swimming toward the top of the screen?
imaligned = 'yes';
%Is the blanking outline of the jelly oriented apex-up?
blankaligned = 'no';
%Is the pressure field oriented correctly (apex-up)?
pfieldaligned = 'no';



%Read apex and manubrium points
%Set Excel file name
apexandmanubrium = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection/apexandmanubrium.xlsx'; %file path and name for Excel file with apex and manubrium points
%Update the ##:## to match cell names in Excel
apexpts = xlsread(apexandmanubrium,'Sheet1','B4:C33');
% apexpts(:,2) = -apexpts(:,2)+(1024*pixscale);
apexpts = apexpts./1000;
manupts = xlsread(apexandmanubrium,'Sheet1','D4:E33');
% manupts(:,2) = -manupts(:,2)+(1024*pixscale);
manupts = manupts./1000;


%Read left and right bell edge points
%Set Excel file name
belledges = 'G:/Jellyfish turns/Aurelia_KNLresults/July 2018/Nov 17_2011_Aurelia 6cm_03_selection/belledges.xlsx'; %file path and name for Excel file with apex and manubrium points
%Update the ##:## to match cell names in Excel
lfpts = xlsread(belledges,'Sheet1','B4:C33');
lfpts = lfpts./1000;
rtpts = xlsread(belledges,'Sheet1','D4:E33');
rtpts = rtpts./1000;



%If you want to inflate the boundary coordinates to try to avoid undefined
%pressure, set snake to yes.  For force calculation directly on the jelly's
%surface, set snake to no.
snake = 'yes';
%Set whether to show the points where force was calculated on the pressure
%contour - yes or no
boundON = 'no';


%Set percentage of bell margin data you want to keep (ex: 0.3 means keep
%outer 30% of bell's bottom). Set to 1 to keep everything
threshold = 0.33;
    
%Set image size in pixels
imwidth = 1104;
imheight = 1024;

%%%%%END MANUAL PARAMETER SETTINGS%%%%%





%Set scales in real units for the image files
xmin = 0;
xmax = xmin+(imwidth*pixscale);
ymin = 0;
ymax = ymin+(imheight*pixscale);

numformat = '%05d';    % format of numbers in file name. '%05d' is five digit format with leading zeros
delimiter = '\t';      % delimiter between columns in blanking file
headerlines = 1;       % number of header lines in velocity file

%Set time step between pressure fields and the number of video frames to be
%analyzed
deltat = increment * 0.001;
numFrames = (last-first)/increment;


%Set a row count index (indicates where the current frame's data will be
%saved in the final data table)
count = 1;
 
%Initialize storage for totals which will be save out later
totpospull = zeros(numFrames,1);
totpospush = zeros(numFrames,1);
totnegpull = zeros(numFrames,1);
totnegpush = zeros(numFrames,1);
totyforce = zeros(numFrames,1);
totxforce = zeros(numFrames,1);
tottorque = zeros(numFrames,1);
tottorquelf = zeros(numFrames,1);
tottorquert = zeros(numFrames,1);


%Initialize storage for x & y forces at each boundary point
forcexall = zeros(numBoundPts,numFrames);
forceyall = zeros(numBoundPts,numFrames);

%Initialize storage for boundary coordinate lists
boundariesx = zeros(numBoundPts,numFrames);
boundariesy = zeros(numBoundPts,numFrames);

%Initialize storage for centroid coordinate lists
surfcentx = zeros(numFrames,1);
surfcenty = zeros(numFrames,1);

%Initialize storage for torque at each boundary point
torqueall = zeros(numBoundPts,numFrames);


%Set up storage for subsetted data
boundxsubset = zeros(numBoundPts,numFrames);
boundysubset = zeros(numBoundPts,numFrames);

forcexsubset = zeros(numBoundPts,numFrames);
forceysubset = zeros(numBoundPts,numFrames);
torquesubset = zeros(numBoundPts,numFrames);

tottorquesubsetlf = zeros(numFrames,1);
tottorquesubsetrt = zeros(numFrames,1);




for filenum = first:increment:last-1    % loop to go through frames and compute forces
    
    close all     % close open figures 
    filenum       % show current filenum in command window
    
    
    
%     if filenum ~= first,
%         prevboundx = boundx;
%         prevboundy = boundy;
%     end

    
    
    
    
    
    apexx = apexpts(count,1);
    apexy = apexpts(count,2);
    manux = manupts(count,1);
    manuy = manupts(count,2);
    
    leftedgex = lfpts(count,1);
    leftedgey = lfpts(count,2);
    rightedgex = rtpts(count,1);
    rightedgey = rtpts(count,2);
    
    
    
    
    
    %%%%%Use these lines to orient image correctly%%%%%
    im = imread([imstub num2str(filenum,numformat) '.' imfileformat]);      % read image file into matrix 'im'.
    if strcmp(imaligned,'no'),
        im = flipud(im);
    end
    %%%%%END image orientation%%%%%
    
    
    
    
    
    press = dlmread([pressstub num2str(filenum,numformat) '.txt'],',',0,0);       %read in pressure field file
    pressXDim = length(find(press(:,2) == press(1,2)));       %find number of x-coordinates of pressure contour (for reshaping field later)
    pressYDim = length(find(press(:,1) == press(1,1)));       %find number of y-coordinates of pressure contour (for reshaping field later)
    if strcmp(pfieldaligned,'no'),
        press(:,2) = ymax/1000 - press(:,2);
    end
    
    
    
    
    
    
    
    
    blank = dlmread([blankstub num2str(filenum,numformat) '.txt'],delimiter,0,0);          % read blanking coordinates into variable 'blank'
    if strcmp(blankaligned, 'no'),
        blank(:,2) = ymax - blank(:,2);
    end
    
    
    %drop adjacent duplicate points - interpolation won't work if two
    %points in a row have the same coordinates.  Here, take difference
    %between a point's coordinate and coordinate of the next point in the
    %list.  If the difference is 0 for both x and y coordinates, drop that
    %point.
    blank = blank(find(abs(blank(:,1)-circshift(blank(:,1),-1)) + abs(blank(:,2)-circshift(blank(:,2),-1)) ~=0),:);
    
    %for interpolation, add the first point to the end of the coordinate
    %list to create a closed loop.
    blank = [blank; blank(1,:)];
    %Interpolate
    blankXtra = interparc(linspace(0,1,numBoundPts+1),blank(:,1),blank(:,2),'pchip');         %Fit a smoothed curve to outline, then pick numBoundPts equi-distant points to be the force calculation boundary
    %Drop the extra point at the end of the list (the repeat for a closed
    %loop)
    blankXtra = blankXtra(1:end-1,:);
    
%     %Turn on below to visually check that the interpolation worked.
%     figure(1)
%     hold on
%     plot(blank(:,1),blank(:,2),'b-')
%     plot(blankXtra(:,1),blankXtra(:,2),'ro')
%     return
%              
    
    blankXtra = blankXtra/1000;     %Convert boundary coordinates from mm to meters
    
    %Convert boundary coordinates to a force calculation boundary
    boundx = blankXtra(:,1);
    boundy = blankXtra(:,2);
    
    
    
    
    
    
    
    %%%%%%%%%%START SNAKE ALGORITHM%%%%%%%%%%%%%%
    if strcmp(snake,'yes'),    
        %Inflate body outline via snake


        %Make binary image of body
        xpix = (1000*boundx)/pixscale;
        ypix = (1000*boundy)/pixscale;
        imBody = poly2mask(xpix,ypix,1024,1024);
%         figure(2)
%         imshow(imBody)
%         axis on
%         hold on
%         %set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)
%         return

        %get distances away from boundary
        Daway = bwdist(imBody);

        imBody = imcomplement(imBody);


        %get distance from boundary
        D = -bwdist(imBody);
        D = D + Daway;
        [rows cols] = size(D);
        [xdim ydim] = meshgrid(1:1:cols,1:1:rows);

        iterations = 8;

        %alpha = elasticity weight, beta = flexibility weight, gamma = step size for moving x and
        %y, and kappa is weight of external gradient
        numpts = size(boundx,1);
        alpha = 0.04;
        beta = 0;%0.2;
        gamma = 0.7;
        kappa = 1.5;

        xs = xpix;
        ys = ypix;

%         return

        for iter = 1:1:iterations,
            %set up elastic force on each point - force keeping points 
            %evenly spread out.  Basically, saying here that we want the
            %average distance between points to be about the same.
            F_elasticx = alpha*2*((circshift(xs,1)-xs)+(circshift(xs,-1)-xs));
            F_elasticy = alpha*2*((circshift(ys,1)-ys)+(circshift(xs,-1)-xs));


            %set up stiffness force on each point - force keeping points 
            %from making too bendy a line.  Only works >2 pts from end.
            %will do this by preventing curvature (1/radcurv) from being
            %too big.  Like above, we'll do this by saying we want the
            %average curvature between adjacent points to not be very 
            %different.  This means we need the curvature at j-1, j, and
            %j+1.

            %curvature 
            m1 = (ys-circshift(ys,1))./(xs-circshift(xs,1));
            m2 = (circshift(ys,-1)-ys)./(circshift(xs,-1)-xs);
            xc = (m1.*m2.*(circshift(ys,1)-circshift(ys,-1))+m2.*(circshift(xs,1)+xs)-m1.*(xs+circshift(xs,-1)))./(2*(m2-m1));
            yc = -(1./m1).*(xc-(circshift(xs,1)+xs)/2)+(circshift(ys,1)+ys)/2;
            rad = sqrt((xs-xc).^2+(ys-yc).^2);

            F_stiff = beta*2*((1./circshift(rad,1)-1./rad)+(1./circshift(rad,-1)-1./rad));


            %Last, set up force to maximize distance from the boundary.
            %If distance at x+1 is bigger, move that way, etc.

            %get distance value at x+1 and x-1
            Dlf = interp2(xdim,ydim,D,xs-2,ys);
            Drt = interp2(xdim,ydim,D,xs+2,ys);
            %get distance value at y+1 and y-1
            Dup = interp2(xdim,ydim,D,xs,ys-2);
            Ddn = interp2(xdim,ydim,D,xs,ys+2);

            %get forces
            F_distx = (kappa/2)*(Drt-Dlf);
            F_disty = (kappa/2)*(Ddn-Dup);



            xs = xs + gamma*(nansum([F_distx,F_stiff,F_elasticx],2));
            ys = ys + gamma*(nansum([F_disty,F_stiff,F_elasticy],2));



%             if iter == 1,
%                 figure(5)
%             end
%             plot(xpix,ypix,'c-','LineWidth',2)
%             hold on
%             plot(xs,ys,'g*')
%             %axis([0 1024 0 1024])
%             axis square
%             hold off
%             pause(0.5)
%             %return

        end
%         return
        boundx = (xs/1000)*pixscale;
        boundy = (ys/1000)*pixscale;
        boundx = double(boundx);
        boundy = double(boundy);
            
            
        %Smooth adjusted outline with an arclength interpolation
        boundx = [boundx; boundx(1)];
        boundy = [boundy; boundy(1)];
        bound = interparc(linspace(0,1,numBoundPts+1),boundx,boundy,'pchip');         %Fit a smoothed curve to outline, then pick numBoundPts equi-distant points to be the force calculation boundary
        bound = bound(1:end-1,:);
        boundx = bound(:,1);
        boundy = bound(:,2);
            
            
%         figure(1)
%         hold on
%         plot(blank(:,1)/1000,blank(:,2)/1000)
%         plot(boundx,boundy,'.')
%         axis equal
%         return


    %%%%%END SNAKE ALGORITHM%%%%%%%%%%%%%
    end
    
    
    





    %Reorder boundary coordinates so apex is first
    
    %Start by finding row number of apex
    
    %Find the distance between each boundary point and the apex point
    disttoapex = ((boundx-apexx).^2+(boundy-apexy).^2).^0.5;
    %Get row number of the point closest to the apex
    [~,apexidx] = min(disttoapex);
    
    %Do the reorder
    boundx = circshift(boundx,-(apexidx-1));
    boundy = circshift(boundy,-(apexidx-1));

    
    
    
    
    
   
    vel = dlmread([vecstub num2str(filenum,numformat) '.txt'],delimiter,headerlines,0);     % read in velocity vector file
    vel(:,1) = vel(:,1)/1000;     % convert x coordinate from mm to meters                          
    vel(:,2) = vel(:,2)/1000;     % convert y coordinate from mm to meters


    
    
    
    
    
    
    %plot image with blanking coordinates overlaid in yellow dots
    figure(1)
    imshow(im,'XData',[xmin/1000 xmax/1000],'YData',[ymin/1000 ymax/1000])
    axis equal
    axis on
    %set(gca,'YDir','normal')
    hold on
    plot(boundx,boundy,'y.')
    

%     return
    






    
    %Plot the pressure contour.
    figure(2)
    contourf(reshape(press(:,1),pressYDim,pressXDim),reshape(press(:,2),pressYDim,pressXDim),reshape(press(:,7),pressYDim,pressXDim),50,'LineColor','none')
    %Set the color palette of the pressure contour.
    if exist('flippedcolors','var'),
        colormap(flippedcolors)
    else
        colormap(jet)
    end
    
    hold on
    fill(blankXtra(:,1),blankXtra(:,2),[0 0 0]);      %Plot jelly silhouette from the blanking boundary data
    if strcmp(boundON,'yes'),
        plot(boundx,boundy,'k.')
    end
    c = colorbar;
    %Set up axes, labels, etc
%     if filenum == first,
%         maxCp = max(abs(press(:,7)));
%     end
%     caxis([-maxCp, maxCp])
    caxis([-3 3])
    xlabel('position [m]')
    ylabel('position [m]')
    c.Label.String = 'Pressure [Pa]';
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);
    axis equal
    
    %axis([(xmin+(xmax-xmin)*0.12)/1000, (xmax-(xmax-xmin)*0.12)/1000, (ymin+(ymax-ymin)*0.12)/1000, (ymax-(ymax-ymin)*0.12)/1000])
    
    set(gca,'YDir','reverse')
%     return
    
    %Save out figure
    if exist('flippedcolors','var'),
        print([savePath '/pressureField_'  num2str(filenum,numformat)], '-dtiffn','-r300');
    else
        print([savePath '/pressureField_jet_'  num2str(filenum,numformat)], '-dtiffn','-r300');
    end
    

  
% pause(2)
% if not(filenum > 20),
%     continue
% end
% return
    
    
    
    
    
    
    

    % calculate pressure at body boundary
    surfpress = griddata(reshape(press(:,1),128,128),reshape(press(:,2),128,128),reshape(press(:,7),128,128),boundx,boundy);

    % define body surface positions that are halfway between the original body
    % boundary points - these surface positions are in the center of the
    % line segments defined by the original boundary points.
    surfposx = (boundx + circshift(boundx,1))./2;
    surfposy = (boundy + circshift(boundy,1))./2;


    surfdx = boundx - circshift(boundx,1);      % compute x-distance between surface points
    surfdy = boundy - circshift(boundy,1);      % compute y-distance between surface points

    darea = sqrt(surfdx.^2 + surfdy.^2);        % compute surface area (per unit depth) between surface points

    surfnormx = surfdy;                       % compute x-component of INWARD pointing vector normal to surface
    surfnormy = -surfdx;                        % compute y-component of INWARD pointing vector normal to surface

    % compute unit vectors
    surfunitnormx = surfnormx./(sqrt(surfnormx.^2 + surfnormy.^2));  
    surfunitnormy = surfnormy./(sqrt(surfnormx.^2 + surfnormy.^2));


    
    
    
    
    
    
    
    %Calculate forces of fluid on body at each boundary point (spacing)*(direction of OUTWARD normal vector)*(pressure) to get units in Pa/m
    forcex(1:size(boundx,1),1) = sqrt(surfdx.^2+surfdy.^2).*-surfunitnormx.*surfpress;
    forcey(1:size(boundy,1),1) = sqrt(surfdx.^2+surfdy.^2).*-surfunitnormy.*surfpress;

    %Get rid of any really outrageously wrong vectors - 10 standard
    %deviations above mean (occur when the gap between boundary points is
    %irregular)
    xfix = find(abs(forcex) > 10*std(abs(forcex(~isnan(forcex))))+mean(abs(forcex(~isnan(forcex)))));
    yfix = find(abs(forcey) > 10*std(abs(forcey(~isnan(forcey))))+mean(abs(forcey(~isnan(forcey)))));
    if ~isempty(xfix),
        forcex(xfix) = NaN;
        forcey(xfix) = NaN;
    end
    if ~isempty(yfix),
        forcex(yfix) = NaN;
        forcey(yfix) = NaN;
    end
    

    
    %Calculate body geometry from blanking boundaries (not force
    %calculation boundaries)
    
    % define body surface positions that are halfway between the original body
    % boundary points - these surface positions are in the center of the
    % line segments defined by the original boundary points.
    surfposbodyx = (blankXtra(:,1) + circshift(blankXtra(:,1),1))./2;
    surfposbodyy = (blankXtra(:,2) + circshift(blankXtra(:,2),1))./2;
    
    
    [geom, iner, cpmo] = polygeom(surfposbodyx, surfposbodyy);                      % compute geometric parameters of object
    
    surfcentx(count,1) = geom(2);                     % compute x-centroid and store in temporal record
    surfcenty(count,1) = geom(3);                     % compute y-centroid and store in temporal record


    %compute moment arm for torque calculation
    armx = surfposx-surfcentx(count,1);   
    army = surfposy-surfcenty(count,1);
    arm = cat(2, armx, army, zeros(size(armx)));                            % moment arm

    surfunitnorm = cat(2, surfunitnormx, surfunitnormy, zeros(size(armx))); % unit normal

    torque = -surfpress.*(cross(arm, surfunitnorm, 2)*[0;0;1]).*darea;     % calculate torque on each surface facet
  

    
    
    %Store boundary coordinates, force, torque at each coordinate
    forcexall(:,count) = forcex;
    forceyall(:,count) = forcey;
    boundariesx(:,count) = boundx;
    boundariesy(:,count) = boundy;
    torqueall(:,count) = torque;
    
    %Store total vertical and horiztonal force
    totyforce(count,1) = nansum(forcey);
    totxforce(count,1) = nansum(forcex);
    
    %Store total torque
    tottorque(count,1) = nansum(torque);
    
    
    
    
    
    %%%%%%%%%%%%ASSIGN LEFT/RIGHT%%%%%%%%%%%%%%
    
    %Find the row number of the manubrium point
    
    %Find the distance between each boundary point and the manubrium point
    disttomanu = ((boundx-manux).^2+(boundy-manuy).^2).^0.5;
    %Get row number of the point closest to the manubrium
    [~,manuidx] = min(disttomanu);
    
    %For ease of computation here, assemble the boundary coordinates into
    %one big matrix
    assembled_bound = [boundx boundy];
    
    %Create a dictionary in a cell structure, where we'll label points
    dictionary = {assembled_bound}; 
    
    %Go through boundary points clockwise.  Label as right until we see
    %the manubrium point, then label left
    for row = 1:1:size(assembled_bound,1),
        if row < manuidx,
            dictionary{1,2}{row,1} = 'right';
        else
            dictionary{1,2}{row,1} = 'left';
        end
    end

    locs(:,count) = dictionary{1,2};

    
    %%%%%%%%%END ASSIGN LEFT/RIGHT%%%%%%%%%%%%
    
    %Sum up torque per left and right side
    tottorquelf(count,1) = nansum(torque(strcmp(locs(:,count),'left')));
    tottorquert(count,1) = nansum(torque(strcmp(locs(:,count),'right')));
    
    
    if not(threshold == 1),
    
        %Keep points on the bottom of bell after a certain threshold length


        %Divide bell boundary into left and right

        %First, need to reorder blankXtra list so that apex is first
        %start by finding distance between each boundary point and apex
        disttoapex = ((blankXtra(:,1)-apexx).^2+(blankXtra(:,2)-apexy).^2).^0.5;
        %Get row number of the point closest to apex
        [~,apexidx] = min(disttoapex);
        %Do the reorder
        blankx = circshift(blankXtra(:,1),-(apexidx-1));
        blanky = circshift(blankXtra(:,2),-(apexidx-1));

        %Next, need to break the lists into left and right
        %start by finding the distance between each boundary point and the manubrium point
        disttomanu = ((blankx-manux).^2+(blanky-manuy).^2).^0.5;
        %Get row number of the point closest to the manubrium
        [~,manuidx] = min(disttomanu);

        %make left and right lists
        rightx = blankx(1:manuidx);
        righty = blanky(1:manuidx);
        leftx = blankx(manuidx+1:end);
        lefty = blanky(manuidx+1:end);

        %Find radius on left and right sides
        rightrad = abs((manuy-apexy)*rightedgex-(manux-apexx)*rightedgey+manux*apexy-manuy*apexx)/sqrt((manuy-apexy)^2+(manux-apexx)^2);
        leftrad = abs((manuy-apexy)*leftedgex-(manux-apexx)*leftedgey+manux*apexy-manuy*apexx)/sqrt((manuy-apexy)^2+(manux-apexx)^2);

        %Get distance between each boundary point and midline
        rightdist = abs((manuy-apexy)*rightx-(manux-apexx)*righty+manux*apexy-manuy*apexx)/sqrt((manuy-apexy)^2+(manux-apexx)^2);
        leftdist = abs((manuy-apexy)*leftx-(manux-apexx)*lefty+manux*apexy-manuy*apexx)/sqrt((manuy-apexy)^2+(manux-apexx)^2);

        %Get cutoff dist for left and right side
        rightcutoff = rightrad-threshold*rightrad;
        leftcutoff = leftrad-threshold*leftrad;

        %Find the points on the bottom of the bell that are LESS THAN this
        %cutoff distance away from midline.  Drop from the subset.

        %Get right and left side force calculation boundaries.
        boundrightx = boundx(strcmp(locs(:,count),'right'));
        boundrighty = boundy(strcmp(locs(:,count),'right'));
        boundleftx = boundx(strcmp(locs(:,count),'left'));
        boundlefty = boundy(strcmp(locs(:,count),'left'));

        %Find distance from midline for each force calculation point
        bounddist = abs((manuy-apexy)*boundx-(manux-apexx)*boundy+manux*apexy-manuy*apexx)/sqrt((manuy-apexy)^2+(manux-apexx)^2);
        bounddistright = bounddist(strcmp(locs(:,count),'right'));
        bounddistleft = bounddist(strcmp(locs(:,count),'left'));

        %Points are listed clockwise, so, for right side, want to keep 
        %everything before the last time dist > cutoff
        rightidx = find(bounddistright > rightcutoff,1,'last');

        %For left side, want to keep everything after the first time dist >
        %cutoff
        leftidx = find(bounddistleft > leftcutoff,1,'first');

        %Get indices for subset calculation boundary
        subsetidx = [rightidx+1:sum(strcmp(locs,'right'))+leftidx-1];


        %     %Turn on the following for visuals useful in debugging
        %     figure(5)
        %     plot(blankXtra(:,1),blankXtra(:,2))
        %     hold on
        % % %     plot(rightx,righty,'r.')
        % % %     plot(leftx,lefty,'g.')
        % % %     plot(rightx(rtmaxidx),righty(rtmaxidx),'bo')
        % % %     plot(leftx(lfmaxidx),lefty(lfmaxidx),'ko')
        % %     plot(boundrightx,boundrighty,'k.')
        % %     plot(boundrightx(1),boundrighty(1),'ko')
        % %     plot(boundrightx(rightidx),boundrighty(rightidx),'ko')
        %     plot(boundx,boundy,'r.')
        %     plot(boundx(subsetidx),boundy(subsetidx),'ko')
        %     set(gca,'YDir','reverse')
        %     return


        boundxsubset(:,count) = boundx;
        boundxsubset(subsetidx,count) = NaN;

        boundysubset(:,count) = boundy;
        boundysubset(subsetidx,count) = NaN;

        forcexsubset(:,count) = forcex;
        forcexsubset(subsetidx, count) = NaN;

        forceysubset(:,count) = forcey;
        forceysubset(subsetidx, count) = NaN;

        torquesubset(:,count) = torque;
        torquesubset(subsetidx, count) = NaN;
        
        tottorquesubsetrt(count,1) = nansum(torque(1:rightidx));
        tottorquesubsetlf(count,1) = nansum(torque(sum(strcmp(locs,'right'))+leftidx:end));

        %%%%%%%%FOR FIGURES ONLY%%%%%%%%%%%
        forcexPlot = forcexsubset(:,count) + circshift(forcexsubset(:,count),-1);
        forcexPlot = forcexPlot(1:2:end);
        forceyPlot = forceysubset(:,count) + circshift(forceysubset(:,count),-1);
        forceyPlot = forceyPlot(1:2:end);
        boundxPlot = boundxsubset(:,count) + circshift(boundxsubset(:,count),-1);
        boundxPlot = boundxPlot/2;
        boundxPlot = boundxPlot(1:2:end);
        boundyPlot = boundysubset(:,count) + circshift(boundysubset(:,count),-1);
        boundyPlot = boundyPlot/2;
        boundyPlot = boundyPlot(1:2:end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %     %Turn on the following for visuals useful in debugging
    %     figure(4)
    %     hold off
    %     imshow(im,'XData',[xmin/1000 xmax/1000],'YData',[ymin/1000 ymax/1000])
    %     hold on
    %     plot(blankXtra(:,1),blankXtra(:,2))
    %     plot(boundx,boundy,'r.')
    %     plot(boundxsubset(:,count),boundysubset(:,count),'ko')
    %     set(gca,'YDir','reverse')
    %     axis equal
    %     pause(1)
    %     
    %     count=count+1;
    %     continue
    end
    
    

    
    
    %Plot total forces acting on boundary coordinates
    figure(3)
%     fill(blankXtra(:,1),blankXtra(:,2),[0 0 0]);      %Plot jelly silhouette from the blanking boundary data
%     hold on
%     %Plot the forces of the fluid acting on the jellyfish at each boundary
%     %coordinate point, scale up each vector by doubling its length
%     if not(threshold==1),
%         quiver(boundxsubset(:,count),boundysubset(:,count),3*forcexsubset(:,count),3*forceysubset(:,count),'AutoScale','off','color',[0.004 0.522 0.443]);
%     else
%         quiver(boundx,boundy,3*forcex,3*forcey,'AutoScale','off','color',[0.004 0.522 0.443]);
%     end
    %%%%%%USE THIS INSTEAD OF ABOVE FOR FIGURE MAKING%%%%%%%
    fill(blankXtra(:,1),blankXtra(:,2),[0.663 0.663 0.663],'LineStyle','none');
    hold on
    if not(threshold==1),
        quiver(boundxPlot,boundyPlot,2*forcexPlot,2*forceyPlot,'AutoScale','off','color',[1 0 0],'LineWidth',2);
    else
        quiver(boundx,boundy,3*forcex,3*forcey,'AutoScale','off','color',[0.004 0.522 0.443]);
    end
    count = count + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.15 0.15 0.7 0.7]);
    axis([min(press(:,1))-0.002 max(press(:,1))+0.002  min(press(:,2))-0.002 max(press(:,2))+0.002])
    daspect([1 1 1])
    
    set(gca,'YDir','reverse')
    
    %Make a scale bar
    
    %Draw a horizontal line of length 0.002 m (= 2*0.001 N/m)(vectors are 
    %plotted on the same scale as the x and y axes.
    line([0.05*range(press(:,1))+min(press(:,1)); 0.05*range(press(:,1))+min(press(:,1))+0.002],[0.1*range(press(:,2))+min(press(:,2)+0.002); 0.1*range(press(:,2))+min(press(:,2))+0.002],'color',[1 0 0],'LineWidth',2);
    %Set dimension label for scale bar 
    str = '1 mN/m';
    %Add dimension label to screen.
    text(0.05*range(press(:,1))+min(press(:,1))+0.0033, 0.1*range(press(:,2))+min(press(:,2))+0.002,str,'FontSize',16);
    xlabel('position [m]')
    ylabel('position [m]')
    


    %Save out figure
    print([savePath '/totalForces_fig_'  num2str(filenum,numformat)], '-dtiffn','-r300');

    continue
    
    
    

    
    

    scale = 3;  % scale factor for plotting Y-FORCE COMPONENT vectors on figure 1
    
    
    indpospull = find(surfunitnormy < 0 & surfpress < 0);   % find surface points where low pressure is pulling animal forward
    indpospush = find(surfunitnormy > 0 & surfpress > 0);   % find surface points where high pressure is pushing animal forward
    indnegpull = find(surfunitnormy > 0 & surfpress < 0);   % find surface points where low pressure is pulling animal backward
    indnegpush = find(surfunitnormy < 0 & surfpress > 0);   % find surface points where high pressure is pushing animal backward



    %Get total forces for all categories above
    if threshold == 1,
        totpospull(count,1) = nansum(forcey(indpospull));
        totpospush(count,1) = nansum(forcey(indpospush));
        totnegpull(count,1) = nansum(forcey(indnegpull));
        totnegpush(count,1) = nansum(forcey(indnegpush));
    else
        totpospull(count,1) = nansum(forceysubset(indpospull,count));
        totpospush(count,1) = nansum(forceysubset(indpospush,count));
        totnegpull(count,1) = nansum(forceysubset(indnegpull,count));
        totnegpush(count,1) = nansum(forceysubset(indnegpush,count));
    end
        
        




    % plot arrow for Y-FORCE COMPONENT according to 4 categories above.
    figure(1)
    if not(threshold == 1),
        hold off
        imshow(im,'XData',[xmin/1000 xmax/1000],'YData',[ymin/1000 ymax/1000])
        axis equal
        axis on
        %set(gca,'YDir','normal')
        hold on
        plot(boundxsubset(:,count),boundysubset(:,count),'y.')
        quiver(boundxsubset(indpospull,count),boundysubset(indpospull,count),zeros(size(forceysubset(indpospull,count))),scale*forceysubset(indpospull,count),0,'Color',[0 1 1],'LineWidth',1);
        quiver(boundxsubset(indpospush,count),boundysubset(indpospush,count),zeros(size(forceysubset(indpospush,count))),scale*forceysubset(indpospush,count),0,'Color',[1 0.5 1],'LineWidth',1);
        quiver(boundxsubset(indnegpull,count),boundysubset(indnegpull,count),zeros(size(forceysubset(indnegpull,count))),scale*forceysubset(indnegpull,count),0,'Color',[0 0 1],'LineWidth',1);
        quiver(boundxsubset(indnegpush,count),boundysubset(indnegpush,count),zeros(size(forceysubset(indnegpush,count))),scale*forceysubset(indnegpush,count),0,'Color',[1 0 0],'LineWidth',1);
    else
        quiver(boundx(indpospull,:),boundy(indpospull,:),zeros(size(forcey(indpospull,:))),scale*forcey(indpospull,:),0,'Color',[0 1 1],'LineWidth',1);
        quiver(boundx(indpospush,:),boundy(indpospush,:),zeros(size(forcey(indpospush,:))),scale*forcey(indpospush,:),0,'Color',[1 0.5 1],'LineWidth',1);
        quiver(boundx(indnegpull,:),boundy(indnegpull,:),zeros(size(forcey(indnegpull,:))),scale*forcey(indnegpull,:),0,'Color',[0 0 1],'LineWidth',1);
        quiver(boundx(indnegpush,:),boundy(indnegpush,:),zeros(size(forcey(indnegpush,:))),scale*forcey(indnegpush,:),0,'Color',[1 0 0],'LineWidth',1);
    end
    
    imMax = 1.024*pixscale;
    axis([xmin, imMax, ymin, imMax])

    
    %Draw a horizontal line of length 0.003 m (= 3*0.001 N/m)(vectors are 
    %plotted on the same scale as the x and y axes.
%     fill([0.03*imMax 0.03*imMax 0.03*imMax+0.022 0.03*imMax+0.022],[0.08*imMax 0.08*imMax+0.01 0.08*imMax+0.01 0.08*imMax],[1 1 1]) 
    line([0.05*imMax; 0.05*imMax+0.003],[0.1*imMax+0.002; 0.1*imMax+0.002],'color','y','LineWidth',2);
    %Set dimension label for scale bar 
    str = '1 mN/m';
    %Add dimension label to screen.
    text(0.05*imMax+0.0033, 0.1*imMax+0.002,str,'FontSize',16,'color','y');
    xlabel('position [m]')
    ylabel('position [m]')
    
    
%     return
    %Save out figure
    print([savePath '/forceMechanism_'  num2str(filenum,numformat)], '-dtiffn','-r300');
        
        
    
    
    
    
    %Increase row count (so next data value is save in the next row)
    count=count+1;
    
    


%     return
end







% %Write data per boundary point to files
dlmwrite([savePath '/boundx.dat'],boundariesx,delimiter);
dlmwrite([savePath '/boundy.dat'],boundariesy,delimiter);
dlmwrite([savePath '/Fx.dat'],forcexall,delimiter);
dlmwrite([savePath '/Fy.dat'],forceyall,delimiter);
dlmwrite([savePath '/torque.dat'],torqueall,delimiter);
dlmwrite([savePath '/torque_RIGHTtotals.dat'],tottorquert,delimiter);
dlmwrite([savePath '/torque_LEFTtotals.dat'],tottorquelf,delimiter);

if not(threshold == 1),
    dlmwrite([savePath '/boundx_subset.dat'],boundxsubset,delimiter);
    dlmwrite([savePath '/boundy_subset.dat'],boundysubset,delimiter);
    dlmwrite([savePath '/Fx_subset.dat'],forcexsubset,delimiter);
    dlmwrite([savePath '/Fy_subset.dat'],forceysubset,delimiter);
    dlmwrite([savePath '/torque_subset.dat'],torquesubset,delimiter);
    dlmwrite([savePath '/torque_subset_RIGHTtotals.dat'],tottorquesubsetrt,delimiter);
    dlmwrite([savePath '/torque_subset_LEFTtotals.dat'],tottorquesubsetlf,delimiter);
end


%Write location markers to file
fileID = fopen([savePath '/locationofpoints.dat'], 'w');
[nrows,ncols]=size(locs);
for row = 1:1:nrows,
    for col = 1:1:ncols,
        if class(locs{row,col}) == 'cell',
            locs{row,col} = locs{row,col}{1,1};
        end
    end
end
for row=1:1:nrows,
    fprintf(fileID,'%s\n',strjoin(locs(row,:),'\t'));
end
fclose(fileID);






if exist('flippedcolors','var'),
    currFr = first;
    v = VideoWriter([savePath '/pressureField.avi'],'Uncompressed AVI');
    open(v);
    for frame = 1:1:numFrames,
        im = imread([savePath '/pressureField_' num2str(currFr,numformat) '.tif']);
        writeVideo(v,im);
        currFr = currFr + increment;
    end
    close(v);
else
    currFr = first;
    v = VideoWriter([savePath '/pressureField_jet.avi'],'Uncompressed AVI');
    open(v);
    for frame = 1:1:numFrames,
        im = imread([savePath '/pressureField_jet_' num2str(currFr,numformat) '.tif']);
        writeVideo(v,im);
        currFr = currFr + increment;
    end
    close(v);
end

% currFr = first;
% v = VideoWriter([savePath '/totalForces.avi'],'Uncompressed AVI');
% open(v);
% for frame = 1:1:numFrames,
%     im = imread([savePath '/totalForces_' num2str(currFr,numformat) '.tif']);
%     writeVideo(v,im);
%     currFr = currFr + increment;
% end
% close(v);
%%%%%%%%%%FIGURE MAKING ONLY%%%%%%%%%%
currFr = first;
v = VideoWriter([savePath '/totalForces_fig.avi'],'Uncompressed AVI');
open(v);
for frame = 1:1:numFrames,
    im = imread([savePath '/totalForces_fig_' num2str(currFr,numformat) '.tif']);
    writeVideo(v,im);
    currFr = currFr + increment;
end
close(v);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% currFr = first;
% v = VideoWriter([savePath '/forceMechanism.avi'],'Uncompressed AVI');
% open(v);
% for frame = 1:1:numFrames,
%     im = imread([savePath '/forceMechanism_' num2str(currFr,numformat) '.tif']);
%     writeVideo(v,im);
%     currFr = currFr + increment;
% end
% close(v);
