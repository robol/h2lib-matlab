it=10:5:20;
h=fopen('risultati.txt', 'w');
fprintf(h,'matrix_dim\t min_block_dim\t iter\t time\t \t res\n');
for n=[500]%,1000,1500,2000,3000,4000,6000,8000]
  for k=ceil(30*1.15.^[0:20])	
       [t,res]=Time(n,k,it);
       for j=1:length(it)
       	fprintf(h,'%d\t \t %d\t \t %d\t %e\t %e\n',n,k,it(j),t(j),res(j));	
       end
    fprintf(h,'\n');
  end
  fprintf(h,'\n');
end
fclose(h);
