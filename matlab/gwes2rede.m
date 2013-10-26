function gwes2rede(gwes_file, rede_file)
% function gwes2rede
%
% Read the results of GWES (http://bioinfo.utu.fi/GWESserver/) from gwes_file
% Convert into the graph format stored in a JSON file
% Write to rede_file

% Read the top pairs
% GWES file has columns:
% SNP1_ID SNP2_ID SNP1_pValue     SNP2_pValue     p-value p-value(-log10) GenoClassNo
% SNP1_Chr        SNP2_Chr        SNP1_Position   SNP2_Position
[id1,id2,uni1,uni2,pval,score,genclass,chrom1,chrom2,pos1,pos2] = textread(gwes_file, ['%s\t%s\t%f\t%f\t%f\%f\t%d\t%d\t%d\t%d\t%d'], 'headerlines','1');
allPairs = load( pairsFile );
[gss, sortIdx] = sort(allPairs.sel.fltGSS, 'descend');
if length(gss) < nPairs
    nPairs = length(gss);
end
topIdx = sortIdx(1:nPairs);

% Output json file
fid = fopen(jsonFile, 'w');
clean_name = regexprep(pairsFile, '\', '/');
fprintf(fid, '{"name": "%s", ', clean_name);

% find top probes
topPairs = allPairs.sel.prb(topIdx,:);
topProbes = sort(unique(reshape(topPairs, 2*nPairs, 1)));
% compute the number of pairs involved with each node
nProbes = length(topProbes);
adjMat = zeros(nProbes, nProbes);
for idx = 1:nPairs
    idx1 = find(topProbes == topPairs(idx,1));
    idx2 = find(topProbes == topPairs(idx,2));
    adjMat(idx1,idx2) = adjMat(idx1,idx2) + 1;
    adjMat(idx2,idx1) = adjMat(idx2,idx1) + 1;
end
degree = sum(adjMat);

% get rs number, chromosome, base pair location for probes
allProbes = load( locFile );
fprintf(fid, '"nodes": [');
for idx = 1:nProbes
    if idx ~= 1
        fprintf(fid, ', ');
    end
    idxProbe = topProbes(idx);
    fprintf(fid, '{"prbCode": "%s", "degree": %1.1f, "prb": %d, ',...
            strtrim(allProbes.probeAnnot.prbCode(idxProbe,:)), degree(idx), idxProbe);
    fprintf(fid, '"rs": "%s", "probe_group": 1, "bp_position": %d, "chrom": %d, "id": %d}',...
            strtrim(allProbes.probeAnnot.prbCode(idxProbe,:)), allProbes.probeAnnot.bp(idxProbe), ...
            allProbes.probeAnnot.chr(idxProbe), idx-1);
end
fprintf(fid, '], ');

% Add GWIS pairs as links
fprintf(fid, '"links": [');
filters = fieldnames(allPairs.sel);
filters(1) = [];
for idx = 1:nPairs
    if idx ~= 1
        fprintf(fid, ', ');
    end    
    idxSource = find(topPairs(idx,1)==topProbes)-1;
    idxTarget = find(topPairs(idx,2)==topProbes)-1;
    idxPair = topIdx(idx);
    fprintf(fid, '{"source": %d, "target": %d, "probe_group": 1, "ct_id": %d, ', ...
            idxSource, idxTarget, idx-1);
    for idxFlt = 1:length(filters)
        if idxFlt ~= 1
            fprintf(fid, ', ');
        end
        fprintf(fid, '"%s": %f', filters{idxFlt}, allPairs.sel.(filters{idxFlt})(idxPair));
    end
    fprintf(fid, '}')
end
fprintf(fid, '], ');


fprintf(fid,'}');
fclose(fid);
