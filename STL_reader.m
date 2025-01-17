%%The program reads an STL file from the folder and outputs a boolean 3D-matrix with points on the inside of the
%%model's boundaries marked as 1, and points outside of the model's boundaries as 0

clc; clear;
stlLines = readlines("sphere.stl"); %inputting the stl file as a
lineCount = length(stlLines); %amount of lines in an stl file
eps = 0.001; %an infinitely small number
res = 20; %resolution of the voxel mesh

%% read STL file
ix = 1; %for-loop iteration counter
for i = 1:lineCount-1
    currentLine = stlLines(i);
    newStr = split(currentLine);
    if(newStr(1) == "vertex")
        stlVortices(ix,1) = str2double(newStr(2));
        stlVortices(ix,2) = str2double(newStr(3));
        stlVortices(ix,3) = str2double(newStr(4));
        ix = ix+1;
    end
end

maxCoords = max(stlVortices); %1x3 row with each axis' maximum given coordinates
minCoords = min(stlVortices); %1x3 row with each axis' minimum given coordinates

polyCount = length(stlVortices)/3; %number of polygons in the stl file

xlinspace = linspace(minCoords(1)-eps,maxCoords(1)+eps,res);
ylinspace = linspace(minCoords(2)-eps,maxCoords(2)+eps,res);
zlinspace = linspace(minCoords(3)-eps,maxCoords(3)+eps,res);
[X,Y,Z] = meshgrid(xlinspace,ylinspace,zlinspace);
%constructs a 3d array with res amount of rows, columns and pages, with elements equally spaced
%between maximum and minimum of each axis
%%
marks = zeros(res,res,res); %constructs a 3d array of zeroes with res amount of pages, rows and columns
parfor k = 1:res
    for j = 1:res
        for i = 1:res
            %ray-polygon intersection counters; max - in the positive direction; min - in the negative direction
            counterxmax = 0;
            counterymax = 0;
            counterzmax = 0;
            counterxmin = 0;
            counterymin = 0;
            counterzmin = 0;
            %point coordinate values to be checked
            checkX = X(i,j,k);
            checkY = Y(i,j,k);
            checkZ = Z(i,j,k);
            for polyCycles = 1:polyCount
                polyAddress = polyCycles*3;
                %poli variables get vertex coordinates for each possible
                %polygon
                polyX1 = stlVortices(polyAddress-2,1);
                polyY1 = stlVortices(polyAddress-2,2);
                polyZ1 = stlVortices(polyAddress-2,3);
                polyX2 = stlVortices(polyAddress-1,1);
                polyY2 = stlVortices(polyAddress-1,2);
                polyZ2 = stlVortices(polyAddress-1,3);
                polyX3 = stlVortices(polyAddress,1);
                polyY3 = stlVortices(polyAddress,2);
                polyZ3 = stlVortices(polyAddress,3);
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(maxCoords(1) + eps,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,maxCoords(1) + eps,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,maxCoords(1) + eps,checkY,checkZ,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,maxCoords(1) + eps,checkY,checkZ,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterxmax = counterxmax+1;
                    end
                end
                
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(checkX,maxCoords(2) + eps,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,checkX,maxCoords(2) + eps,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,checkX,maxCoords(2) + eps,checkZ,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,checkX,maxCoords(2) + eps,checkZ,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterymax = counterymax+1;
                    end
                end
                
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(checkX,checkY,maxCoords(3) + eps,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,maxCoords(3) + eps,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,maxCoords(3) + eps,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,maxCoords(3) + eps,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterzmax = counterzmax+1;
                    end
                end
                
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(checkX,checkY,minCoords(1) - eps,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,minCoords(1) - eps,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,minCoords(1) - eps,checkY,checkZ,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,minCoords(1) - eps,checkY,checkZ,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterxmin = counterxmin+1;
                    end
                end
                
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(checkX,minCoords(2) - eps,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,checkX,minCoords(2) - eps,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,checkX,minCoords(2) - eps,checkZ,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,checkX,minCoords(2) - eps,checkZ,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterymin = counterymin+1;
                    end
                end
                
                SVolume1 = SignedVolume(checkX,checkY,checkZ,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                SVolume2 = SignedVolume(checkX,checkY,minCoords(3) - eps,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                if((SVolume1*SVolume2 <= 0))
                    SVolume3 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,minCoords(3) - eps,polyX1,polyY1,polyZ1,polyX2,polyY2,polyZ2);
                    SVolume4 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,minCoords(3) - eps,polyX2,polyY2,polyZ2,polyX3,polyY3,polyZ3);
                    SVolume5 = SignedVolume(checkX,checkY,checkZ,checkX,checkY,minCoords(3) - eps,polyX3,polyY3,polyZ3,polyX1,polyY1,polyZ1);
                    if (((sign(SVolume3) == sign(SVolume4))||((SVolume3*SVolume4) == 0))&&((sign(SVolume4) == sign(SVolume5))||((SVolume3*SVolume4) == 0))&&((sign(SVolume5) == sign(SVolume3))||((SVolume5*SVolume3) == 0)))
                        counterzmin = counterzmin+1;
                    end
                end
            end
            if ((isOdd(counterxmax)==1)&&(isOdd(counterymax)==1)&&(isOdd(counterzmax)==1))
                marks(i,j,k) = 1;
                continue
            end
            if ((isOdd(counterxmin)==1)&&(isOdd(counterymin)==1)&&(isOdd(counterzmin)==1))
                marks(i,j,k) = 1;
            end
        end
    end
end
%% this part of the program gives a visual representation of the voxelization
f1 = figure;
f2 = figure;
figure(f1);
fv = isosurface(X,Y,Z,marks);
p = patch(fv);
isonormals(X,Y,Z,marks,p)
p.FaceColor = 'blue';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3);
camlight
lighting gouraud

%% HCE
Tinit = 300.0; % [K]
T = marks .* Tinit;
% T = T + Tinit;

time = 0.0;
tau = 0.01;
tmax = 10.0;
h = zlinspace(2) - zlinspace(1);

ambT = 300.0;
preheatTime = 0.2; % [s]
dwellTime = 0.1; % [s]
layer = 1.5*h; % [m]
currentLayer = 0.0; % [m]
currentTime = preheatTime + dwellTime + 1; % [s]
tempMarks = marks * 0.0;
maxZMarks = 2;
gs = 100;

while (time < tmax)
    if ( currentTime > (preheatTime + dwellTime) )
        currentLayer = currentLayer + layer;
        currentTime = 0.0;
        
        for ix = 1:res
            for iy = 1:res
                for iz = 1:res
                    if (Z(1,1,iz) < currentLayer)
                        tempMarks(ix,iy,iz) = marks(ix,iy,iz);
                    end
                    
                end
            end
        end
        
        for ix = 1:res
            for iy = 1:res
                for iz = 1:res
                    if (tempMarks(ix,iy,iz) == 1 && iz >= maxZMarks)
                        maxZMarks = iz;
                    end
                    
                end
            end
        end
        
    end
    if (currentTime < preheatTime)
        T(:,:,maxZMarks) = 2000 .* tempMarks(:,:,maxZMarks); % [K]
        Tprev = T;
    else
        T(:,:,maxZMarks) = tempMarks(:,:,maxZMarks).*(T(:,:,maxZMarks-1)+1.5*h*ambT)/(1+2*h);
        
    end
    
    for GS = 1:gs
        for iz = 2:maxZMarks-1
            thermalDiffusivity = (-10.906*T(floor(res/2),floor(res/2),iz) + 17707.8)/900; %Сталь 20
            for iy = 2:res-1
                for ix = 2:res-1
                    Tsum = T(ix-1,iy,iz) + T(ix+1,iy,iz) + T(ix,iy-1,iz) + T(ix,iy+1,iz) + T(ix,iy,iz+1) + T(ix,iy,iz-1);
                    temp = thermalDiffusivity .* tau .* marks(ix,iy,iz) ./ h.^2;
                    T(ix,iy,iz) = marks(ix,iy,iz)*( Tsum * temp + Tprev(ix,iy,iz) ) ./ (1 + temp .* 6);
                end
            end
        end
    end
    figure(f2);
    xslice = mean([maxCoords(1);minCoords(1)]);
    yslice = mean([maxCoords(2);minCoords(2)]);
    zslice = [];
    slice(X,Y,Z,T,xslice,yslice,zslice)
    %     view([90 0])
    title(['Time is ',num2str(time),' s']);
    %     colormap jet
    newmap = jet;                    %starting map
    ncol = size(newmap,1);           %how big is it?
    zpos = 1;                        %only the 1st position
    newmap(zpos,:) = [1 1 1];        %set that position to white
    colormap(newmap);                %activate it
    colorbar
    caxis([100 1500]);
    %     pause(0.1);
    drawnow;   
    
    Tprev = T;
    time = time + tau;
    currentTime = currentTime + tau;
end

%% functions
function SV = SignedVolume(cx,cy,cz,p1x,p1y,p1z,p2x,p2y,p2z,p3x,p3y,p3z)
%the function returns a signed volume of a matrix input
checkCollision = [p1x-cx, p1y-cy, p1z-cz; p2x-cx, p2y-cy, p2z-cz; p3x-cx, p3y-cy, p3z-cz];
SV = det(checkCollision)/6;
end

function tf = isOdd(x) %determines whether an input is odd with an output of 1 if input is odd and 0 if the input is even
tf = mod(x,2) == 1;
end