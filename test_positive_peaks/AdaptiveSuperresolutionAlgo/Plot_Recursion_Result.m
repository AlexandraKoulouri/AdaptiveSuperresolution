function Plot_Recursion_Result(Hnew,x_new,y_new,x,y,Ind_nz_cell,xi,yi) 


figure; set(gcf, 'Units','centimeters', 'Position',[5 5 9 7])

i=1;

for i= 1:length(x_new)

  if ~isempty(Hnew{i,1}) 
     for jj = 1:size(Hnew{i,1},1)
                hold on
                patch(x_new{i}(Hnew{i,1}(jj,:)),y_new{i}(Hnew{i,1}(jj,:)),[227 227 227]./255,'EdgeColor',[160 160 160]./255);
        end

     end
end     

%       hold on
hold on        
plot(xi,yi,'xr','markersize',8)
      
 i = 1;
 while i<=length(x)      
       
     plot(x{i}(Ind_nz_cell{i}),y{i}(Ind_nz_cell{i}),'ob','markersize',2)%,'linewidth',1.3);
     i = i+1;
 end
    %axis([0 1 0 1]) 
   % title(['Iteration m = ',num2str(iter-1)],'FontSize',8); 
axis square
box on 
axis xy 
%axis([0.45 0.5 0.273 0.295])
set(gcf, 'renderer', 'painter');
set(gcf,'PaperUnits','centimeters') 


