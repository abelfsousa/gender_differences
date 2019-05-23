#!/bin/bash

#convert network from cytoscape to gml format


mean_fpkm=$1
edge_list=$2
file_name=$3

cat $mean_fpkm $edge_list |
awk 'BEGIN{FS="\t";id=1;}{
if(NF==2){geneexp[$1]=$2;};
if(NF>2){
    if( ($0 !~/weight/) && ($5!="" && $6!="") ){
        key=$5" -- "$6;
        if($3>0){edgeList[key] = $3;
            if(!($5 in hashgid)){gid[id]=$5;hashgid[$5]=id;id=id+1;}
            if(!($6 in hashgid)){gid[id]=$6;hashgid[$6]=id;id=id+1;}
        }
    }
}}END{
    #print nodes
    print "graph"
    print "["
    print "directed 0"
    for(g=1;g<=length(gid);g++){
        print "\tnode"
        print "\t["
        print "\tid "g
        print "\tlabel \""gid[g]"\""
        gname = gid[g]
        if((geneexp[gname]<1)||(!(gname in geneexp))){print "\tgroup 1";print "\tfill \"#fcae91\""}
        if((geneexp[gname]>1)&&(geneexp[gname]<10)){print "\tgroup 2";print "\tfill \"#de2d26\""}
        if(geneexp[gname]>10){print "\tgroup 3";print "\tfill \"#a50f15\""}
        print "\tborder \"#aba7c4\""
        print "\t]"
    }
    for(k in edgeList){
        split(k, arr, " -- ")
        print "\tedge"
        print "\t["
        print "\tsource "hashgid[arr[1]]
        print "\ttarget "hashgid[arr[2]]
        cr = edgeList[k] * 10;
        print "\tvalue "cr;
        print "\t]"
    }
    print "]"
}' > $file_name