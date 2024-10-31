function bd=coco_bd_table(run)
bdc=coco_bd_read(run);
bd=cell2table(bdc(2:end,1:end-1),'VariableNames',bdc(1,1:end-1));
end
