#include "xcorr_io.h"
#include "xcorr_postprocess.h"
#include <csignal>
extern "C" {
#include "xcorr2_args.h"
}
int main(int argc,char **argv) {
    std::signal(SIGXFSZ,SIG_IGN);
    try {
        bool geo=argc>1 && std::string(argv[1])=="geo";
        std::vector<std::string> names={"freq_xcorr.dat"};
        if(geo)names.insert(names.end(),{"azi_offset.grd","rng_offset.grd","azi_offset_ll.grd","rng_offset_ll.grd"});
        OutputTransaction outputs(names);
        CheckedFile file(outputs.path("freq_xcorr.dat"));
        if(getenv("UNBUFFERED"))setvbuf(file.get(),nullptr,_IONBF,0);
        for(int y=100;y<=600;y+=100)for(int x=100;x<=600;x+=100)
            if(std::fprintf(file.get(),"%d 1 %d 2 80 10\n",x,y)<0)throw io_error("write test correlation");
        file.close();
        st_xcorr xc{};xc.do_geocode=geo;xc.do_blockmedian=getenv("RAW_GRID")==nullptr;
        xc.xsearch=xc.ysearch=128;xc.ri=2;xc.snr_thr=getenv("EMPTY_FILTER")?100:10;xc.psnr_thr=5;
        xc.nxl=xc.nyl=6;
        std::vector<int> axis={100,200,300,400,500,600};
        run_postprocess(xc,"master.PRM",axis,axis,outputs);
        outputs.publish();return 0;
    } catch(const std::exception &e) {std::fprintf(stderr,"xcorr_cc: %s\n",e.what());return 1;}
}
