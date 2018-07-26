function Clusters = MergeEdges(Edges,Magnitude,Direction)
    %Constants to export sometime
    thetaThr = 100;
    magThr = 1200;
    
    ValidIds = unique([Edges(:,2) ; Edges(:,3)]);
    EncRing = [ValidIds,(1:size(ValidIds,1))'];
    
    %convert Edges to local mapping
    Loops = length(ValidIds);
    parfor k = 2:3
        ThisCol = Edges(:,k);
        for j = 1:Loops
            ThisCol(ThisCol(:) == ValidIds(j)) = j;
        end
        Edges(:,k) = ThisCol
    end

%     for j = 1:length(ValidIds)
%         Edges(ValidIds(j) == Edges(:,2),2) = j;
%         Edges(ValidIds(j) == Edges(:,3),3) = j; 
%     end
    
    %Reshape the Magnitude and Directions of the arrays
    tmin = ArraytoList(Direction);
    tmin = tmin(ValidIds);
    tmax = tmin;
    mmin = ArraytoList(Magnitude);
    mmin = mmin(ValidIds);
    mmax = mmin;
    
    %Create the unionfind vector which is pre allocated for speed
    SimpleUF = [(1:length(ValidIds))',ones(length(ValidIds),1)];
    
    for i = 1:size(Edges,1)
        ida = Edges(i,2);
        idb = Edges(i,3);
        
        ida = IgetRepresentative(SimpleUF,ida); %gets rep
        idb = IgetRepresentative(SimpleUF,idb); %gets rep
        
        if(ida == idb) %It's already connected!
            continue;
        end
        
        sza = SimpleUF(ida,2); %Get the size of tree a
        szb = SimpleUF(idb,2); %Get the size of tree b
        
        tmina = tmin(ida); tmaxa = tmax(ida); %finds the max and mins of a
        tminb = tmin(idb); tmaxb = tmax(idb); %finds the max and mins of b
        
        costa = tmaxa-tmina; %Intermediate cost value
        costb = tmaxb-tminb; %Intermediate cost value
        
        %Makes sure that the angles aren't more than +/- PI
        bshift = mod2pi((tminb+tmaxb)/2,(tmina+tmaxa)/2)-(tminb+tmaxb)/2;
        
        tminab = min(tmina, tminb+bshift); %Theta min
        tmaxab = max(tmaxa, tmaxb+bshift); %Theta max
        
        if(tmaxab-tminab > 2*pi)
            tmaxab = tminab + 2*pi;
        end
        
        mminab = min(mmin(ida), mmin(idb)); %Mag min
        mmaxab = max(mmax(ida), mmax(idb)); %Mag max
        
        costab = (tmaxab - tminab); %Intermediate cost value for a and b
        
        %Magic Values that I need to understand more :)
        Value1 = (costab <= min(costa,costb) + (thetaThr/(sza+szb)));
        Value2 = (mmaxab-mminab) <= min(mmax(ida)-mmin(ida),...
            mmax(idb)-mmin(idb)) + (magThr/(sza+szb));
        
        if(Value1 && Value2)
            [SimpleUF, idab] = IconnectNodes(SimpleUF, ida, idb);
            
            tmin(idab) = tminab; %Sets the minimum theta
            tmax(idab) = tmaxab; %Sets the maximum theta
            
            mmin(idab) = mminab; %Sets the minimum mag
            mmin(idab) = mmaxab; %Sets the maximum mag
        end
    end
    
    %Export the clusters
    Clusters = ExportClusters(SimpleUF,Magnitude, Edges,EncRing);
end

% function localIds = Cvt2LocMany(Mapping,EdgeMatrix)
% localIds(Mapping(:,1) == EdgeMatrix(:,2),2) = Mapping(:,2); 
% end

function localId = Cvt2Loc(Mapping,EdgeId)
    localId = Mapping(EdgeId == Mapping(:,1),2);
end

function globalId = Cvt2Pic(Mapping,localId)
    globalId = Mapping(localId == Mapping(:,2),1);
end

% Gets the representative of the node
function root = IgetRepresentative(UFArray,NodeId)
    if(UFArray(NodeId,1) == NodeId) %If it is it's own rep return
        root = NodeId;              %No changes
    else
        root = UFArray(NodeId,1);
    end
end

% Gets the representative of the node
function [root,UpdatedArray] = getRepresentative(UFArray,NodeId)
    if(UFArray(NodeId,1) == NodeId) %If it is it's own rep return
        root = NodeId;              %No changes
    else
        root = getRepresentative(UFArray,UFArray(NodeId,1)); %Recurse
        UFArray(NodeId,1) = root; %Flatten the tree
    end
UpdatedArray = UFArray;   %Return the updated array
end

%connects and merges the two trees together
function [UFArray,root] = IconnectNodes(UFArray, aId,bId)

    aRoot = IgetRepresentative(UFArray,aId); %Get rep of a
    bRoot = IgetRepresentative(UFArray,bId); %Get rep of b

    if(aRoot==bRoot) %It's already connected!
        root=aRoot;  %Return the root
        return;
    end
    
    if(UFArray(aRoot,2) > UFArray(bRoot,2)) %Larger tree wins!
        %Add the sizes together
        UFArray(aRoot,2) = UFArray(aRoot,2) + UFArray(bRoot,2);
        
        logicArr = UFArray(UFArray == bRoot);
        for j = 1:length(UFArray)
            if(logicArr)
                UFArray(j,1) = aRoot;
            end
        end        
        
        root=aRoot; %Return the new root
        return;
    else
        %Add the sizes together
        UFArray(bRoot,2) = UFArray(aRoot,2) + UFArray(bRoot,2);
        
        logicArr = UFArray(UFArray == aRoot);
        for j = 1:length(UFArray)
            if(logicArr)
                UFArray(j,1) = bRoot;
            end
        end
        
        
        root=bRoot; %Return the new root
        return;
    end
end

%connects and merges the two trees together
function [NewUFArray,root] = connectNodes(UFArray, aId,bId)

    [aRoot, UFArray] = getRepresentative(UFArray,aId); %Get rep of a
    [bRoot, UFArray] = getRepresentative(UFArray,bId); %Get rep of b
    NewUFArray = UFArray; %copy the array

    if(aRoot==bRoot) %It's already connected!
        root=aRoot;  %Return the root
        return;
    end
    
    if(UFArray(aRoot,2) > UFArray(bRoot,2)) %Larger tree wins!
        NewUFArray(bRoot,1) = aRoot; %Set the new root
        
        %Add the sizes together
        NewUFArray(aRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2);
        
        root=aRoot; %Return the new root
        return;
     else
        NewUFArray(aRoot,1) = bRoot; %Set the new root
        
        %Add the sizes together
        NewUFArray(bRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2);
        
        root=bRoot; %Return the new root
        return;
    end
end

function longArray = ArraytoList(Array)
Width = size(Array,2);
Height  = size(Array,1);

longArray = zeros(1,Width*Height);
for i = 1:Height
    StartIdx = ((i-1) * Width)+1;
    EndIdx   = (StartIdx + Width)-1;
    longArray(1,StartIdx:EndIdx) = Array(i,:);
end
end

%Formatting the clusters as a list with points to make it easier later
function ClusterList = ExportClusters(UF_Array,Magnitude,Edges,Mapping)
    %Need to export these constants
    MinCluster = 4;

    %find clusters that have more than the MinSeg
    Valid_Clusters = UF_Array((UF_Array(:,2) >= MinCluster),1);
    

    %Create a logical array for faster indexing / display
    logical_arr = ismember(UF_Array(:,1),Valid_Clusters);
    
    ClusterList = zeros(size(Edges,1),4);  %Empty matrix for clusters
    
    for i = 1:size(Edges,1)-1 %loops through all the edges
        EdgeCluster = UF_Array(Edges(i,3),1); %Gets cluster #
        EdgeMag = Magnitude(Edges(i,3));      %Gets magnitude
        EdgeX = Edges(i,4);                   %Gets X coord
        EdgeY = Edges(i,5);                   %Gets Y coord
        ClusterList(i,:) = [EdgeX,EdgeY,EdgeMag,EdgeCluster];
    end

    ClusterList = sortrows(ClusterList,4);
end