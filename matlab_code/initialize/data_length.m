function [IY1, IM1, ID1, IH1, IY2, IM2, ID2, IH2] = data_length(IYB, IMB, IDB, IHB, IYE, IME, IDE, IHE)

    start_date = datetime(IYB, IMB, IDB, IHB, 0, 0);
    end_date = datetime(IYE, IME, IDE, IHE, 0, 0);
    time_interval = hours(1);
    
    time_sequence = (start_date:time_interval:end_date)';
    IY1 = year(time_sequence);
    IM1 = month(time_sequence);
    ID1 = day(time_sequence);
    IH1 = hour(time_sequence);
    

    start_date = datetime(IYB, IMB, IDB);
    end_date = datetime(IYE, IME, IDE);
    dates = start_date:end_date;

    IY2 = year(dates)';
    IM2 = month(dates)';
    ID2 = day(dates)';
    IH2 = zeros(length(IY2), 1);
end