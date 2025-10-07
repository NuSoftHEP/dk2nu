#include "calcLocationWeights.h"
#include "dk2nu.h"
#include "dkmeta.h"

#include "pybind11/pybind11.h"

#include "pybind11/native_enum.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TVector3.h"

#include <memory>
#include <tuple>

namespace py = pybind11;

using namespace bsim;

namespace pybind11 {
namespace detail {

// adapted from
// https://pybind11.readthedocs.io/en/latest/advanced/cast/custom.html
template <> struct type_caster<TVector3> {

  PYBIND11_TYPE_CASTER(TVector3, io_name("Sequence[float]",
                                         "tuple[float, float, float]"));

  bool load(handle src, bool /*convert*/) {
    // Check if handle is a Sequence
    if (!py::isinstance<py::sequence>(src)) {
      return false;
    }
    auto seq = py::reinterpret_borrow<py::sequence>(src);
    // Check if exactly two values are in the Sequence
    if (seq.size() != 3) {
      return false;
    }
    // Check if each element is either a float or an int
    for (auto item : seq) {
      if (!py::isinstance<py::float_>(item) &&
          !py::isinstance<py::int_>(item)) {
        return false;
      }
    }
    value.SetXYZ(seq[0].cast<double>(), seq[1].cast<double>(),
                 seq[2].cast<double>());
    return true;
  }
};

} // namespace detail
} // namespace pybind11

class dk2nuFileReader {

  std::unique_ptr<TChain> dk2nu_ch;
  std::unique_ptr<TTreeReader> dk2nu_rdr;
  std::unique_ptr<TTreeReaderValue<Dk2Nu>> dk2nu;

  std::unique_ptr<TChain> dkmeta_ch;
  std::unique_ptr<TTreeReader> dkmeta_rdr;
  std::unique_ptr<TTreeReaderValue<DkMeta>> dkmeta;

  bool isdone;

public:
  dk2nuFileReader(std::vector<std::string> const &infiles) {

    dk2nu_ch = std::make_unique<TChain>("dk2nuTree");
    dkmeta_ch = std::make_unique<TChain>("dkmetaTree");

    for (auto const &f : infiles) {
      dk2nu_ch->Add(f.c_str());
      dkmeta_ch->Add(f.c_str());
    }

    dk2nu_rdr = std::make_unique<TTreeReader>(dk2nu_ch.get());
    dk2nu = std::make_unique<TTreeReaderValue<Dk2Nu>>(*dk2nu_rdr, "dk2nu");
    dkmeta_rdr = std::make_unique<TTreeReader>(dkmeta_ch.get());
    dkmeta = std::make_unique<TTreeReaderValue<DkMeta>>(*dkmeta_rdr, "dkmeta");
  }

  Dk2Nu const *_first() {
    dk2nu_rdr->Restart();
    isdone = false;
    return _next();
  }

  Dk2Nu const *_next() {
    if (dk2nu_rdr->Next()) {
      return &**dk2nu;
    } else {
      return nullptr;
    }
  }

  pybind11::object first() {
    dk2nu_rdr->Restart();
    isdone = false;
    return next();
  }

  pybind11::object next() {
    if (dk2nu_rdr->Next()) {
      return py::cast(**dk2nu);
    } else {
      return py::none();
    }
  }

  pybind11::object first_meta() {
    dkmeta_rdr->Restart();
    isdone = false;
    return next_meta();
  }
  pybind11::object next_meta() {
    if (dkmeta_rdr->Next()) {
      return py::cast(**dkmeta);
    } else {
      return py::none();
    }
  }
};

class dk2nuFileReader_sentinel {};

class dk2nuFileReader_iter {
  std::reference_wrapper<dk2nuFileReader> dkfr;
  pybind11::object curr_event;

public:
  dk2nuFileReader_iter(dk2nuFileReader &_dkfr) : dkfr(_dkfr) {
    curr_event = dkfr.get().first();
  }
  void operator++() { curr_event = dkfr.get().next(); }
  pybind11::object const &operator*() { return curr_event; }
  bool operator!=(dk2nuFileReader_sentinel const &) const {
    return !curr_event.is(py::none());
  }
  bool operator==(dk2nuFileReader_sentinel const &) const {
    return curr_event.is(py::none());
  }
};

dk2nuFileReader_iter begin(dk2nuFileReader &dkfr) {
  return dk2nuFileReader_iter(dkfr);
}

class dk2nuFileReader_metaiter {
  std::reference_wrapper<dk2nuFileReader> dkfr;
  pybind11::object curr_entry;

public:
  dk2nuFileReader_metaiter(dk2nuFileReader &_dkfr) : dkfr(_dkfr) {
    curr_entry = dkfr.get().first_meta();
  }
  void operator++() { curr_entry = dkfr.get().next_meta(); }
  pybind11::object const &operator*() { return curr_entry; }
  bool operator!=(dk2nuFileReader_sentinel const &) const {
    return !curr_entry.is(py::none());
  }
  bool operator==(dk2nuFileReader_sentinel const &) const {
    return curr_entry.is(py::none());
  }
};

dk2nuFileReader_metaiter begin_meta(dk2nuFileReader &dkfr) {
  return dk2nuFileReader_metaiter(dkfr);
}

dk2nuFileReader_sentinel end(dk2nuFileReader &) {
  return dk2nuFileReader_sentinel();
}

PYBIND11_MODULE(pydk2nu, m) {
  py::module::import("ROOT");

  py::native_enum<flgbitval>(m, "flgbitval", "enum.Enum")
      .value("kFlgOverflow", kFlgOverflow)
      .value("kMaskReserved", kMaskReserved)
      .value("kMaskUser", kMaskUser)
      .export_values()
      .finalize();

  py::native_enum<dkproc>(m, "dkproc", "enum.Enum")
      .value("dkp_unknown", dkp_unknown)
      .value("dkp_k0l_nuepimep", dkp_k0l_nuepimep)
      .value("dkp_k0l_nuebpipem", dkp_k0l_nuebpipem)
      .value("dkp_k0l_numupimmup", dkp_k0l_numupimmup)
      .value("dkp_k0l_numubpipmum", dkp_k0l_numubpipmum)
      .value("dkp_kp_numumup", dkp_kp_numumup)
      .value("dkp_kp_nuepi0ep", dkp_kp_nuepi0ep)
      .value("dkp_kp_numupi0mup", dkp_kp_numupi0mup)
      .value("dkp_kp_numubmum", dkp_kp_numubmum)
      .value("dkp_kp_nuebpi0em", dkp_kp_nuebpi0em)
      .value("dkp_kp_numubpi0mum", dkp_kp_numubpi0mum)
      .value("dkp_mup_nusep", dkp_mup_nusep)
      .value("dkp_mum_nusep", dkp_mum_nusep)
      .value("dk_pip_numumup", dk_pip_numumup)
      .value("dk_pim_numubmum", dk_pim_numubmum)
      .value("dk_mum_capture", dk_mum_capture)
      .value("dkp_maximum", dkp_maximum)
      .value("dkp_other", dkp_other)
      .export_values()
      .finalize();

  py::class_<NuRay>(m, "NuRay")
      .def(py::init<double, double, double, double, double>())
      .def_readwrite("px", &NuRay::px)
      .def_readwrite("py", &NuRay::py)
      .def_readwrite("pz", &NuRay::pz)
      .def_readwrite("E", &NuRay::E)
      .def_readwrite("wgt", &NuRay::wgt)
      .def("__str__", &NuRay::AsString, py::arg("opt") = "");

  py::class_<Decay>(m, "Decay")
      .def(py::init<>())
      .def_readwrite("norig", &Decay::norig)
      .def_readwrite("ndecay", &Decay::ndecay)
      .def_readwrite("ntype", &Decay::ntype)
      .def_readwrite("vx", &Decay::vx)
      .def_readwrite("vy", &Decay::vy)
      .def_readwrite("vz", &Decay::vz)
      .def_readwrite("pdpx", &Decay::pdpx)
      .def_readwrite("pdpy", &Decay::pdpy)
      .def_readwrite("pdpz", &Decay::pdpz)
      .def_readwrite("ppdxdz", &Decay::ppdxdz)
      .def_readwrite("ppdydz", &Decay::ppdydz)
      .def_readwrite("pppz", &Decay::pppz)
      .def_readwrite("ppenergy", &Decay::ppenergy)
      .def_readwrite("ppmedium", &Decay::ppmedium)
      .def_readwrite("ptype", &Decay::ptype)
      .def_readwrite("muparpx", &Decay::muparpx)
      .def_readwrite("muparpy", &Decay::muparpy)
      .def_readwrite("muparpz", &Decay::muparpz)
      .def_readwrite("mupare", &Decay::mupare)
      .def_readwrite("necm", &Decay::necm)
      .def_readwrite("nimpwt", &Decay::nimpwt)
      .def_readwrite("sumnimpwt2", &Decay::sumnimpwt2)
      .def("__str__", &Decay::AsString, py::arg("opt") = "");

  py::class_<Ancestor>(m, "Ancestor")
      .def(py::init<>())
      .def_readwrite("pdg", &Ancestor::pdg)
      .def_readwrite("startx", &Ancestor::startx)
      .def_readwrite("starty", &Ancestor::starty)
      .def_readwrite("startz", &Ancestor::startz)
      .def_readwrite("startt", &Ancestor::startt)
      .def_readwrite("startpx", &Ancestor::startpx)
      .def_readwrite("startpy", &Ancestor::startpy)
      .def_readwrite("startpz", &Ancestor::startpz)
      .def_readwrite("stoppx", &Ancestor::stoppx)
      .def_readwrite("stoppy", &Ancestor::stoppy)
      .def_readwrite("stoppz", &Ancestor::stoppz)
      .def_readwrite("polx", &Ancestor::polx)
      .def_readwrite("poly", &Ancestor::poly)
      .def_readwrite("polz", &Ancestor::polz)
      .def_readwrite("nucleus", &Ancestor::nucleus)
      .def_readwrite("parIndex", &Ancestor::parIndex)
      .def_readwrite("proc", &Ancestor::proc)
      .def_readwrite("ivol", &Ancestor::ivol)
      .def_readwrite("imat", &Ancestor::imat)
      .def("__str__", &Ancestor::AsString, py::arg("opt") = "");

  py::class_<TgtExit>(m, "TgtExit")
      .def(py::init<>())
      .def_readwrite("tvx", &TgtExit::tvx)
      .def_readwrite("tvy", &TgtExit::tvy)
      .def_readwrite("tvz", &TgtExit::tvz)
      .def_readwrite("tpx", &TgtExit::tpx)
      .def_readwrite("tpy", &TgtExit::tpy)
      .def_readwrite("tpz", &TgtExit::tpz)
      .def_readwrite("tptype", &TgtExit::tptype)
      .def_readwrite("tgen", &TgtExit::tgen)
      .def("__str__", &TgtExit::AsString, py::arg("opt") = "");

  py::class_<Traj>(m, "Traj")
      .def(py::init<>())
      .def_readwrite("trkx", &Traj::trkx)
      .def_readwrite("trky", &Traj::trky)
      .def_readwrite("trkz", &Traj::trkz)
      .def_readwrite("trkpx", &Traj::trkpx)
      .def_readwrite("trkpy", &Traj::trkpy)
      .def_readwrite("trkpz", &Traj::trkpz)
      .def("__str__", &Traj::AsString, py::arg("opt") = "");

  auto decay_through_point = [](Dk2Nu const &dk2nu,
                                TVector3 const &point_in_beam_frame_cm) {
    double enu_GeV, wght_m2;
    int ec = calcEnuWgt(&dk2nu, point_in_beam_frame_cm, enu_GeV, wght_m2);
    double wght_cm2 = wght_m2 * 1E-4;
    if (ec) {
      std::stringstream ss;
      ss << "Failed to calculate decay energy and weight. dk2nu error "
            "code: "
         << ec;
      ss << "\nDk2Nu object:\n" << dk2nu.AsString();
      ss << "\nDecay point: [" << point_in_beam_frame_cm[0] << ", "
         << point_in_beam_frame_cm[1] << ", " << point_in_beam_frame_cm[2]
         << "] cm." << std::endl;
      throw std::runtime_error(ss.str());
    }
    return std::make_tuple(enu_GeV, wght_cm2 * dk2nu.decay.nimpwt / M_PI);
  };

  py::class_<Dk2Nu>(m, "Dk2Nu")
      .def(py::init<>())
      .def_readwrite("job", &Dk2Nu::job)
      .def_readwrite("potnum", &Dk2Nu::potnum)
      .def_readwrite("jobindx", &Dk2Nu::jobindx)
      .def_readwrite("decay", &Dk2Nu::decay)
      .def_readwrite("nuray", &Dk2Nu::nuray)
      .def_readwrite("ancestor", &Dk2Nu::ancestor)
      .def_readwrite("ppvx", &Dk2Nu::ppvx)
      .def_readwrite("ppvy", &Dk2Nu::ppvy)
      .def_readwrite("ppvz", &Dk2Nu::ppvz)
      .def_readwrite("tgtexit", &Dk2Nu::tgtexit)
      .def_readwrite("traj", &Dk2Nu::traj)
      .def_readwrite("flagbits", &Dk2Nu::flagbits)
      .def_readwrite("vint", &Dk2Nu::vint)
      .def_readwrite("vdbl", &Dk2Nu::vdbl)
      .def("__str__", &Dk2Nu::AsString, py::arg("opt") = "")
      .def("decay_through_point", decay_through_point,
           py::arg("point_in_beam_frame_cm"), R"(
Returns the (neutrino Energy (GeV), total flux weight in /cm^2])
  for this decay to send a neutrino through a 1 m diameter circular window at 
  point_in_beam_frame_cm.

N.B. The importance weight from the beam/target simulation (Decay::nimpwt) 
  is included in the returned flux weight, so the returned weight is the only
  weight you need to make flux predictions.
)");

  py::class_<Location>(m, "Location")
      .def(py::init<>())
      .def_readwrite("x", &Location::x)
      .def_readwrite("y", &Location::y)
      .def_readwrite("z", &Location::z)
      .def_readwrite("name", &Location::name)
      .def("__str__", &Location::AsString, py::arg("opt") = "");

  py::class_<DkMeta>(m, "DkMeta")
      .def(py::init<>())
      .def_readwrite("job", &DkMeta::job)
      .def_readwrite("pots", &DkMeta::pots)
      .def_readwrite("beamsim", &DkMeta::beamsim)
      .def_readwrite("physics", &DkMeta::physics)
      .def_readwrite("physcuts", &DkMeta::physcuts)
      .def_readwrite("tgtcfg", &DkMeta::tgtcfg)
      .def_readwrite("horncfg", &DkMeta::horncfg)
      .def_readwrite("dkvolcfg", &DkMeta::dkvolcfg)
      .def_readwrite("beam0x", &DkMeta::beam0x)
      .def_readwrite("beam0y", &DkMeta::beam0y)
      .def_readwrite("beam0z", &DkMeta::beam0z)
      .def_readwrite("beamhwidth", &DkMeta::beamhwidth)
      .def_readwrite("beamvwidth", &DkMeta::beamvwidth)
      .def_readwrite("beamdxdz", &DkMeta::beamdxdz)
      .def_readwrite("beamdydz", &DkMeta::beamdydz)
      .def_readwrite("location", &DkMeta::location)
      .def_readwrite("vintnames", &DkMeta::vintnames)
      .def_readwrite("vdblnames", &DkMeta::vdblnames)
      .def("__str__", &DkMeta::AsString, py::arg("opt") = "");

  py::class_<dk2nuFileReader>(m, "dk2nuFileReader")
      .def(py::init<std::vector<std::string> const &>())
      .def("first", &dk2nuFileReader::first)
      .def("next", &dk2nuFileReader::next)
      .def(
          "decays",
          [](dk2nuFileReader &s) {
            return py::make_iterator(begin(s), end(s));
          },
          py::keep_alive<0, 1>())
      .def("first_meta", &dk2nuFileReader::first_meta)
      .def("next_meta", &dk2nuFileReader::next_meta)
      .def(
          "metas",
          [](dk2nuFileReader &s) {
            return py::make_iterator(begin_meta(s), end(s));
          },
          py::keep_alive<0, 1>())
      .def("decay_all_through_point",
           [=](dk2nuFileReader &rdr, TVector3 const &point_in_beam_frame_cm) {
             std::vector<std::tuple<int, double, double>> data;

             auto dk = rdr._first();

             while (dk) {

               auto enu_wght = decay_through_point(*dk, point_in_beam_frame_cm);
               data.emplace_back(dk->decay.ntype, std::get<0>(enu_wght),
                                 std::get<1>(enu_wght));

               dk = rdr._next();
             }

             py::array_t<double> a;
             a.resize(std::array<size_t, 2>{data.size(), 3});

             // Obtain mutable access to the array
             auto r = a.mutable_unchecked<2>();

             for (size_t i = 0; i < data.size(); i++) {
               r(i, 0) = std::get<0>(data[i]);
               r(i, 0) = std::get<1>(data[i]);
               r(i, 0) = std::get<2>(data[i]);
             }

             return a;
           },
           py::arg("point_in_beam_frame_cm"), R"(
Returns the (neutrino species, neutrino Energy (GeV), total flux weight in /cm^2])
  for each decay in the file to send a neutrino through a 1 m diameter circular 
  window at point_in_beam_frame_cm and returns the result as a single Nx4 numpy
  array.

N.B. The importance weight from the beam/target simulation (Decay::nimpwt) 
  is included in the returned flux weight, so the returned weight is the only
  weight you need to make flux predictions.
)");
}