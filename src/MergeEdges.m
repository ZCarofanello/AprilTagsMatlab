function Clusters = MergeEdges(Edges,Magnitude,Direction)
    %Constants to export sometime
    thetaThr = 200;
    magThr = 1200;
    
    %Get the width and height of the iamge
    width = size(Magnitude,2);
    height = size(Magnitude,1);
    
    %Reshape the Magnitude and Directions of the arrays
    tmin = reshape(Direction,1,[]);
    tmax = reshape(Direction,1,[]);
    mmin = reshape(Magnitude,1,[]);
    mmax = reshape(Magnitude,1,[]);
    
    %Create the unionfind vector which is pre allocated for speed
    SimpleUF = [(1:width*height)', ones(1,width*height)'];
    
    for i = 1:size(Edges,1)
        ida = Edges(i,2);
        idb = Edges(i,3);
        
        [ida, SimpleUF] = getRepresentative(SimpleUF,ida); %gets rep
        [idb, SimpleUF] = getRepresentative(SimpleUF,idb); %gets rep
        
        if(ida == idb) %It's already connected!
            continue;
        end
        
        sza = SimpleUF(ida,2); %Get the size of tree a
        szb = SimpleUF(idb,2); %Get the size of tree b
        
        tmina = tmin(ida); tmaxa = tmax(ida); %finds the max and mins of a
        tminb = tmin(idb); tmaxb = tmax(idb); %finds the max and mins of b
        
        costa = tmaxa-tmina; %Intermediate cost value
        costb = tmaxb-tminb; %Intermediate cost value
        
        bshift = mod2pi((tminb+tmaxb)/2, (tmina+tmaxa)/2) - (tminb+tmaxb)/2; %Makes sure that the angles aren't more than +/- PI
        
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
        Value2 = (mmaxab-mminab) <= min(mmax(ida)-mmin(ida), mmax(idb)-mmin(idb)) + (magThr/(sza+szb));
        
        if(Value1 && Value2)
            [SimpleUF, idab] = connectNodes(SimpleUF, ida, idb);
            
            tmin(idab) = tminab;
            tmax(idab) = tmaxab;
            
            mmin(idab) = mminab;
            mmin(idab) = mmaxab;
        end
    end
    %Export the clusters
    Clusters = ExportClusters(SimpleUF,Magnitude, Edges);
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
        NewUFArray = FlattenTree(bRoot,aRoot,NewUFArray);
        NewUFArray(aRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2); %Add the sizes together
        root=aRoot; %Return the new root
        return;
     else
        NewUFArray(aRoot,1) = bRoot; %Set the new root
        NewUFArray = FlattenTree(aRoot,bRoot,NewUFArray);
        NewUFArray(bRoot,2) = NewUFArray(aRoot,2) + NewUFArray(bRoot,2); %Add the sizes together
        root=bRoot; %Return the new root
        return;
    end
end

function FixedTree = FlattenTree(OldRoot,NewRoot,UF_Array)
Children_Idx = FindChildren(UF_Array,OldRoot);
UF_Array(Children_Idx,1) = NewRoot;
FixedTree = UF_Array;
end

function Nodes = FindChildren(UF_Array, ParentIdx)
Nodes = find(UF_Array(:,1) == ParentIdx);
end

%Formatting the clusters as a list with points to make it easier later
function ClusterList = ExportClusters(UF_Array,Magnitude,Edges)
    %Need to export these constants
    MinCluster = 4;

    %find clusters that have more than the MinSeg
    Valid_Clusters = UF_Array((UF_Array(:,2) >= MinCluster),1);
    
    for k = 1:length(Valid_Clusters)
        root = getRepresentative(UF_Array,Valid_Clusters(k));
        UF_Array = FlattenTree(Valid_Clusters(k),root,UF_Array);
    end

    %Create a logical array for faster indexing / display
    logical_arr = ismember(UF_Array(:,1),Valid_Clusters);
    
    ClusterList = [];
    FirstEntry = true;
    for i = 1:size(Edges,1)-1
        if(logical_arr(Edges(i,2)))
            EdgeCluster = UF_Array(Edges(i,2),1);
            EdgeMag = Magnitude(Edges(i,5),Edges(i,4));
            EdgeX = Edges(i,5);
            EdgeY = Edges(i,4);
            if(FirstEntry)
                ClusterList = [EdgeX,EdgeY,EdgeMag,EdgeCluster];
                FirstEntry = false;
            else
                ClusterList = [ClusterList;EdgeX,EdgeY,EdgeMag,EdgeCluster];
            end
        end
    end
    Cluster_Num = unique(ClusterList(:,4)); %Gets each unique cluster

    ClusterList = sortrows(ClusterList,4);
end