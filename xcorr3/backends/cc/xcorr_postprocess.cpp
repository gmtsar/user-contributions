#include "xcorr_postprocess.h"
#include "xcorr_io.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/wait.h>
#include <fcntl.h>
extern "C" {
#include "xcorr2.h"
#include "xcorr2_args.h"
}

/* No shell pipeline: every command has its own checked status. stderr remains
 * visible, and all fixed-name GMTSAR scratch files stay in a private directory. */
static void command(const std::vector<std::string> &args, const std::string &cwd,
                    const std::string &stdout_path = "") {
    std::vector<char *> argv;
    for (const auto &arg : args) argv.push_back(const_cast<char *>(arg.c_str()));
    argv.push_back(nullptr);
    int fd = -1;
    if (!stdout_path.empty()) {
        fd = open(stdout_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
        if (fd < 0) throw io_error("open command output " + stdout_path);
    }
    pid_t child = fork();
    if (child < 0) { int saved=errno; if(fd>=0)close(fd); errno=saved; throw io_error("fork " + args[0]); }
    if (child == 0) {
        if (chdir(cwd.c_str()) != 0) _exit(126);
        if (fd >= 0) {
            if (dup2(fd, STDOUT_FILENO) < 0) _exit(126);
            close(fd);
        }
        execvp(argv[0], argv.data());
        _exit(127);
    }
    if (fd >= 0) close(fd);
    int status; pid_t waited;
    do { waited=waitpid(child,&status,0); } while(waited<0 && errno==EINTR);
    if (waited<0) throw io_error("wait for " + args[0]);
    std::string label=args[0]+(args.size()>1 ? " "+args[1] : "");
    if (!WIFEXITED(status) || WEXITSTATUS(status)!=0)
        throw std::runtime_error("geocoding stage failed: " + label +
            (WIFEXITED(status) ? " (exit "+std::to_string(WEXITSTATUS(status))+")" : " (terminated by signal)"));
}

struct PrmOwner {
    prm_handler handler{};
    explicit PrmOwner(const char *path) {
        if (!prm_open(&handler,path)) throw std::runtime_error(handler.error);
    }
    ~PrmOwner() { prm_close(&handler); }
    double get(const char *key) {
        double value;
        if (!prm_get_f64(&handler,key,&value)) throw std::runtime_error(handler.error);
        return value;
    }
};

static void verify_grid(const std::string &grid, OutputTransaction &outputs) {
    OutputTransaction::require_product(outputs.path(grid));
    const std::string info=outputs.path(grid+".info");
    command({"gmt","grdinfo",grid,"-Cn","-L0"},outputs.directory(),info);
    std::ifstream input(info);
    double xmin,xmax,ymin,ymax,zmin,zmax;
    if (!(input>>xmin>>xmax>>ymin>>ymax>>zmin>>zmax) ||
        !std::isfinite(xmin)||!std::isfinite(xmax)||!std::isfinite(ymin)||!std::isfinite(ymax)||
        !std::isfinite(zmin)||!std::isfinite(zmax) || xmin>=xmax || ymin>=ymax || zmin>zmax)
        throw std::runtime_error("invalid or all-NaN geocoding grid: " + grid);
}

void run_postprocess(const st_xcorr &xc, const char *master_prm,
                     const std::vector<int> &xpos, const std::vector<int> &ypos,
                     OutputTransaction &outputs) {
    if (!xc.do_geocode) return;
    char *resolved=realpath("trans.dat",nullptr);
    if (!resolved) throw io_error("-geocode requires readable trans.dat");
    std::string trans=resolved; std::free(resolved);
    OutputTransaction::require_product(trans);
    if (access(trans.c_str(),R_OK)!=0) throw io_error("read trans.dat");
    if (symlink(trans.c_str(),outputs.path("trans.dat").c_str())!=0) throw io_error("stage trans.dat");
    // Preserve proj_ra2ll.csh's existing optional gauss_* sampling hints.
    DIR *hints=opendir(".");
    if (!hints) throw io_error("read geocoding sampling hints");
    while (dirent *entry=readdir(hints)) {
        if (std::strncmp(entry->d_name,"gauss_",6)!=0) continue;
        char *full=realpath(entry->d_name,nullptr);
        if (!full) { closedir(hints); throw io_error("resolve geocoding sampling hint"); }
        int result=symlink(full,outputs.path(entry->d_name).c_str());
        int saved=errno; std::free(full);
        if (result!=0) {closedir(hints);errno=saved;throw io_error("stage geocoding sampling hint");}
    }
    closedir(hints);
    if (xpos.size()<2 || ypos.size()<2)
        throw std::runtime_error("-geocode requires at least two sample positions on each axis");
    PrmOwner prm(master_prm);
    double prf=prm.get("PRF"), velocity=prm.get("SC_vel"), radius=prm.get("earth_radius"),
           height=prm.get("SC_height"), rate=prm.get("rng_samp_rate");
    if (prf<=0 || velocity<=0 || radius<=0 || rate<=0 || 1.0+height/radius<=0)
        throw std::runtime_error("invalid PRM pixel-size parameters for -geocode");
    double azi_size=velocity/std::sqrt(1.0+height/radius)/prf, rng_size=299792458.0/rate/2.0;
    if (!std::isfinite(azi_size) || !std::isfinite(rng_size) || azi_size<=0 || rng_size<=0)
        throw std::runtime_error("nonfinite or nonpositive pixel size for -geocode");
    double max_rng=static_cast<double>(xc.xsearch)/xc.ri-2.0, max_azi=xc.ysearch-2.0;
    std::ifstream fin(outputs.path("freq_xcorr.dat"));
    if (!fin) throw std::runtime_error("cannot read staged correlation output");
    CheckedFile fa(outputs.path("pot_azi.xyz")), fr(outputs.path("pot_rng.xyz"));
    std::string line; size_t total=0, kept=0;
    while (std::getline(fin,line)) {
        double xp,xo,yp,yo,cc,ps; std::istringstream row(line); std::string extra;
        if (!(row>>xp>>xo>>yp>>yo>>cc>>ps) || (row>>extra) ||
            !std::isfinite(xp)||!std::isfinite(xo)||!std::isfinite(yp)||!std::isfinite(yo)||!std::isfinite(cc)||!std::isfinite(ps))
            throw std::runtime_error("malformed staged correlation row " + std::to_string(total+1));
        ++total;
        if (cc>xc.snr_thr && ps>xc.psnr_thr && xo>-max_rng && xo<max_rng && yo>-max_azi && yo<max_azi) {
            if (!std::isfinite(yo*azi_size) || !std::isfinite(xo*rng_size))
                throw std::runtime_error("pixel-to-metre conversion overflow");
            // Keep the existing filtering, range sign and six-significant-digit conversion.
            if (std::fprintf(fa.get(),"%.0f %.0f %.6g\n",xp,yp,yo*azi_size)<0 ||
                std::fprintf(fr.get(),"%.0f %.0f %.6g\n",xp,yp,-xo*rng_size)<0)
                throw io_error("write geocoding coordinates");
            ++kept;
        }
    }
    if (fin.bad() || !fin.eof()) throw std::runtime_error("read staged correlation output failed");
    fa.close(); fr.close();
    if (!kept) throw std::runtime_error("-geocode: no points passed correlation/peak-SNR/search-window filters");
    char region[128],increment[128];
    std::snprintf(region,sizeof(region),"-R%d/%d/%d/%d",xpos.front(),xpos.back(),ypos.front(),ypos.back());
    std::snprintf(increment,sizeof(increment),"-I%.10g/%.10g",double(xpos.back()-xpos.front())/(xpos.size()-1),double(ypos.back()-ypos.front())/(ypos.size()-1));
    for (const std::string comp : {"azi","rng"}) {
        std::string xyz="pot_"+comp+".xyz", grid=comp+"_offset.grd";
        if (xc.do_blockmedian) {
            std::string median="pot_"+comp+"_median.xyz",raw=comp+"_raw.grd";
            command({"gmt","blockmedian",xyz,region,increment},outputs.directory(),outputs.path(median));
            OutputTransaction::require_product(outputs.path(median));
            command({"gmt","xyz2grd",median,region,increment,"-G"+raw},outputs.directory());
            OutputTransaction::require_product(outputs.path(raw));
            command({"gmt","grdfilter",raw,"-Fm3","-Dp","-Nr","-G"+grid},outputs.directory());
        } else command({"gmt","xyz2grd",xyz,region,increment,"-G"+grid},outputs.directory());
        verify_grid(grid,outputs);
        command({"proj_ra2ll.csh","trans.dat",grid,comp+"_offset_ll.grd"},outputs.directory());
        verify_grid(comp+"_offset_ll.grd",outputs);
    }
    std::fprintf(stderr,"xcorr_cc: geocoding stages completed: %zu/%zu points retained; awaiting output publication\n",kept,total);
}
