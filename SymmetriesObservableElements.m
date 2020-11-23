clear all 

% Info for non-matlab users like me
%
% cell: A cell array is a data type with indexed data containers called cells, 
% unique(A): only unique values of A -> remove duplicates
% A == B: returns logical array
% find(X): returns vector of inidces of each non-zero element of X


load('Square_Lattice_Reduced_Embeddings_Order10.mat');

order=4;

double_touch=false;

% gnum is array of graph keys
% consturct empty cell elements of dim 1
relevant_links=cell(size(gnum,1),1);
%       where each cell can contain any type of data.
% size(A, dim): Returns dim-th entry of size vector.
% numel(A): Number of array elements
ground_states=false(size(gnum,1),1);
partyhard_info=cell(size(gnum,1),1);

% iterate over indices of all graph keys
% for i0=1:size(Graphs_upto_Order12,1)
for i0=1:size(gnum,1)



% I have no idea what SmallGrphs is.
% We accees graph g by index i0
% G=Graphs_upto_Order12{i0};
G=SmallGraphs{i0};


% only consider graphs to a given order
if numedges(G)<order+1 

% create block-cut tree graph
% each component of a tree a biconnected component and cut vertex of G.
% ind is a vector of nodes mapping G nodes to tree nodes. 
% Node i of graph G is mapped to index ind(i) of tree
% Only ind is used!
[tree,ind]=bctree(G);

% numel(unique(ind)) = number of biconnected components
A=zeros(numel(unique(ind)),1);
B=unique(ind);

% iterate over indices of biconneted components
for i=1:numel(unique(ind))
    if numel(find(B(i)==ind))==1
        A(i,1)=find(B(i)==ind);
    end
end

A(A==0)=[];

pos_double_Edges=zeros(size(A,1)*(size(A,1)-1)/2,1);

Edges=table2array(G.Edges);
Edges=sort(Edges,2);

m=0;

for i=1:size(A,1)
    for j=1:i-1
        hilfe_1=sort([A(i,1),A(j,1)]);
        m=m+1;
        if sum(sum(hilfe_1==Edges,2)==2)>0
            pos_double_Edges(m,1)=find(sum(hilfe_1==Edges,2)==2);
        end
    end
end

pos_double_Edges(pos_double_Edges==0)=[];

if isempty(pos_double_Edges)
    pos_double_Edges=[];
end

nodes_double_Edges=zeros(size(pos_double_Edges,1),numnodes(G));

for i=1:size(pos_double_Edges,1)
    G1=G;
    G1=rmedge(G1,pos_double_Edges(i,1));
    bins = conncomp(G1);
    if numel(unique(bins))>1
        nodes_double_Edges(i,bins==1)=0;
        nodes_double_Edges(i,bins==2)=1;
    end
end

nodes_double_Edges(sum(nodes_double_Edges,2)==0,:)=[];

minorder=[];
for k1=1:numnodes(G)
    for k2=k1:numnodes(G)
        minorder=[minorder;k1,k2,false];
    end
end
    


for i=1:size(minorder,1)
    m=0;
    m=m+numedges(G)-size(nodes_double_Edges,1);
    for j=1:size(nodes_double_Edges,1)
        if nodes_double_Edges(j,minorder(i,1))==nodes_double_Edges(j,minorder(i,2)) && double_touch
            m=m+2;
        else
            m=m+1;
        end
    end
    minorder(i,3)=m<order+1;
end
test=minorder;
        
minorder=minorder(logical(minorder(:,3)),:);


if ~isempty(minorder)
% relevant_links_order12{i0,1}=minorder;
relevant_links{i0,1}=minorder;
if mod(i0,1000)==0
i0
end
end


n=(numedges(G)-size(nodes_double_Edges,1))+2*size(nodes_double_Edges,1);
% ground_states_order12(i0,1)=n<order+1;
ground_states(i0,1)=n<order+1;

    
end

end



% Reduce number of hopping elements by using symmetries, here for all
% graphs with 11 and 12 edges

% for i=1:size(Graphs_upto_Order12,1)
for i=1:size(SmallGraphs,1)
    partyhard_info_=zeros(numnodes(SmallGraphs{i,1}));

%         if ~isempty(relevant_links_order12{i})
        if ~isempty(relevant_links{i})
%             Links=relevant_links_order12{i};
            Links=relevant_links{i};
            Links=Links(:,1:2);
            
% G=Graphs_upto_Order12{i};
G=SmallGraphs{i};


node_distances=distances(G);

numel_nearest_nodes=zeros(numnodes(G),numedges(G));
for h1=1:numnodes(G)
    for h2=1:numedges(G)
        numel_nearest_nodes(h1,h2)=numel(nearest(G,h1,h2));
    end
end

degree_nodes=degree(G);

perm_1_2=[2,1];

for j1=1:size(Links,1)
    if j1<size(Links,1)+1
G1=G;

partyhard_info_(Links(j1,1),Links(j1,2))=sub2ind([size(partyhard_info_,1),size(partyhard_info_,2)],Links(j1,1),Links(j1,2));
partyhard_info_(Links(j1,2),Links(j1,1))=sub2ind([size(partyhard_info_,1),size(partyhard_info_,2)],Links(j1,1),Links(j1,2));

nodename=zeros(numnodes(G),1);
nodename(Links(j1,1),1)=1;
nodename(Links(j1,2),1)=1;
G1.Nodes.Labels=nodename;
bool=true(size(Links,1),1);
degree_1=sort([degree_nodes(Links(j1,1)),degree_nodes(Links(j1,2))]);
[min_degree_1,I_min_degree_1]=min([degree_nodes(Links(j1,1)),degree_nodes(Links(j1,2))]);
neighbours_min_degree_1=sort(degree_nodes(neighbors(G,Links(j1,I_min_degree_1))));
neighbours_max_degree_1=sort(degree_nodes(neighbors(G,Links(j1,perm_1_2(I_min_degree_1)))));
neighbours_degree_1=sort([neighbours_min_degree_1.',neighbours_max_degree_1.']);


G1_nearest=zeros(numedges(G),2);
for h1=1:numedges(G)
G1_nearest(h1,:)=sort([numel_nearest_nodes(Links(j1,1),1),numel_nearest_nodes(Links(j1,2),1)]);
end


for j2=j1+1:size(Links,1)
    degree_2=sort([degree_nodes(Links(j2,1)),degree_nodes(Links(j2,2))]);
    G2_nearest=zeros(numedges(G),2);
for h1=1:numedges(G)
G2_nearest(h1,:)=sort([numel_nearest_nodes(Links(j2,1),1),numel_nearest_nodes(Links(j2,2),1)]);
end
if isequal(degree_1,degree_2) && node_distances(Links(j1,1),Links(j1,2))==node_distances(Links(j2,1),Links(j2,2)) && isequal(G1_nearest,G2_nearest)
    [min_degree_2,I_min_degree_2]=min([degree_nodes(Links(j2,1)),degree_nodes(Links(j2,2))]);
neighbours_min_degree_2=sort(degree_nodes(neighbors(G,Links(j2,I_min_degree_2))));
neighbours_max_degree_2=sort(degree_nodes(neighbors(G,Links(j2,perm_1_2(I_min_degree_2)))));
neighbours_degree_2=sort([neighbours_min_degree_2.',neighbours_max_degree_2.']);


if isequal(neighbours_degree_1,neighbours_degree_2)
    
G2=G;


nodename=zeros(numnodes(G),1);
nodename(Links(j2,1),1)=1;
nodename(Links(j2,2),1)=1;
G2.Nodes.Labels=nodename;

bool(j2,1)=~isisomorphic(G1,G2,'NodeVariables','Labels');

if not(bool(j2,1))
partyhard_info_(Links(j2,1),Links(j2,2))=sub2ind([size(partyhard_info_,1),size(partyhard_info_,2)],Links(j1,1),Links(j1,2));
partyhard_info_(Links(j2,2),Links(j2,1))=sub2ind([size(partyhard_info_,1),size(partyhard_info_,2)],Links(j1,1),Links(j1,2));
end

end

end

end

Links=Links(bool,:);

    end
end

% relevant_links_order12{i}=Links;
relevant_links{i}=Links;
if mod(i,1000)==0
i
end
        end
        partyhard_info{i,1}=partyhard_info_;
end

% Patrick's inefficient stuff

out_str = "#graph_number  hopping_from hopping_to hopping_symmetry_number global_symmetry_number";
out_str = sprintf(out_str + '\n');

partyhard_info{2}
partyhard_info{1}
glob_sym_numb = 0;
for idx=1:size(partyhard_info)
    % just consider upper triangular matrix, as this symmetry is already considered in my other porgram
    sym_matrix = partyhard_info{idx};
    
    for n=1:numel(sym_matrix) % N^2 elements
       % determine symmetry number
       sym_numb = numel(find(sym_matrix == n));
       if sym_numb ~= 0
           %increment global symmetry number
           glob_sym_numb = glob_sym_numb + sym_numb;
           % extract subscripts from linear index
           [from, to] = ind2sub(size(sym_matrix), n);
           %[from, to] = [from, to]-1;
           % create line for hopping of graph with symmetry number
           % subscript shift
           graphname = "_" + gnum{idx} + "_" + "order" + num2str(numedges(SmallGraphs{idx})) + ".cfg";
           out_str = out_str + graphname + ' ' + num2str(from-1) + ' ' + num2str(to-1) + ' ' + num2str(sym_numb) + ' ' + num2str(glob_sym_numb);
           out_str = sprintf(out_str + '\n');
           if from ~= to
                glob_sym_numb = glob_sym_numb + sym_numb;
                out_str = out_str + graphname + ' ' + num2str(to-1) + ' ' + num2str(from-1) + ' ' + num2str(sym_numb) + ' ' + num2str(glob_sym_numb);
                out_str = sprintf(out_str + '\n');
           end
       end  
    end
end

% write out
filename = "~/Repos/GraphSymmetrien/output/_0_symObsGraphsList_order" + num2str(order) + ".cfg";
fid = fopen(filename,'wt');
fprintf(fid, out_str);
fclose(fid);
