

% write P and Plabels into a user specified delimited file

function write_tfa_matrix(tfas,tfa_ids,exp_labels,filename)

fid = fopen(filename,'w');
if length(tfa_ids) ~= size(tfas,1),
    warndlg('The length of tfa_ids must match the number of rows in the TFA matrix. Could not write TFA to file')
    return;
end

% print header

if length(exp_labels)>0,
  fprintf(fid,'Name');
  for i=1:length(exp_labels),
    fprintf(fid,'\t%s',char(exp_labels(i)));
  end
  fprintf(fid,'\n');
end
for i=1:length(tfa_ids),
    fprintf(fid,'%s',char(tfa_ids(i)));
    for j=1:size(tfas,2),
        fprintf(fid,'\t%d',tfas(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);

        