function savepdf(fig,fname)
    set(fig,'Units','inches')
    scrp = get(fig,'Position');
    set(fig,'PaperPosition',[0 0 scrp(3:4)],'PaperSize',[scrp(3:4)]);
    print(fname,'-dpdf','-vector')
end