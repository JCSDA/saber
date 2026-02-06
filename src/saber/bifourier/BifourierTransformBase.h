/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once


#include <map>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include <boost/noncopyable.hpp>

#include "eckit/mpi/Comm.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransformParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(BifourierTransformParameters, Parameters)

 public:
  // FFT backend
  oops::Parameter<std::string> fftBackend{"fft backend", "fftw", this};

  // Truncation type
  oops::Parameter<std::string> truncationType{"truncation type", "arome", this};

  // Truncation factor (2.0 for a linear grid)
  oops::Parameter<double> truncationFactor{"truncation factor", 2.0, this};

  // Sub-ellipses half-width for calibration (AROME default is 1.5)
  oops::Parameter<double> dwGlb{"sub-ellipses half-width", -1.0, this};

  // Skip tests
  oops::Parameter<bool> skipTests{"skip tests", false, this};

  // Spectral tests tolerance
  oops::Parameter<double> specTolerance{"spectral tolerance", 1.0e-9, this};
};

// -----------------------------------------------------------------------------

class BifourierTransformBase : public util::Printable,
                               private boost::noncopyable {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransformBase";}

  // Constructor
  BifourierTransformBase(const oops::GeometryData &,
                         const std::string &,
                         const oops::Variables &,
                         const BifourierTransformParameters &);

  // Destructor
  virtual ~BifourierTransformBase()
    {}

  // Forward FFT
  virtual void gp2sp(std::vector<double> &,
                     atlas::util::Metadata &) const = 0;

  // Inverse FFT
  virtual void sp2gp(std::vector<double> &,
                     const atlas::util::Metadata &) const = 0;

  // Forward FFT adjoint
  virtual void gp2spAdj(std::vector<double> &,
                        const atlas::util::Metadata &) const = 0;

  // Inverse FFT adjoint
  virtual void sp2gpAdj(std::vector<double> &,
                        atlas::util::Metadata &) const = 0;

  // Forward FFT normalization
  virtual double gp2spNorm(const size_t &) const
    {return 1.0;}

  // Inverse FFT normalization
  virtual double sp2gpNorm(const size_t &) const
    {return 1.0;}

  // Forward FFT adjoint normalization
  virtual double gp2spAdjNorm(const size_t &) const
    {return 1.0;}

  // Inverse FFT adjoint normalization
  virtual double sp2gpAdjNorm(const size_t &) const
    {return 1.0;}

  // Communication vectors
  virtual const size_t & recvSize() const = 0;
  virtual const std::vector<int> & recvCounts() const = 0;
  virtual const std::vector<int> & recvDispls() const = 0;

  // Non-virtual methods

  // Accessors

  // Geometry data
  const oops::GeometryData & geometryData() const
    {return gdata_;}

  // UIDs
  const std::string & gridUid() const
    {return gridUid_;}
  const std::string & specUid() const
    {return specUid_;}

  // Grid size
  const size_t & nx() const
    {return nx_;}
  const size_t & ny() const
    {return ny_;}

  // Grid cell sizes
  const double & dx() const
    {return dx_;}
  const double & dy() const
    {return dy_;}

  // Grid mean latitude
  const double & meanLat() const
    {return meanLat_;}

  // Spectral field size
  const size_t & ns() const
    {return ns_;}

  // Spectral field global size
  const size_t & nsGlb() const
    {return nsGlb_;}

  // Spectral field size on each task
  const std::vector<size_t> & nsPerTask() const
    {return nsPerTask_;}

  // Communication vectors
  const std::vector<int> & sCounts() const
    {return sCounts_;}
  const std::vector<int> & sDispls() const
    {return sDispls_;}
  const std::vector<size_t> & sMapping() const
    {return sMapping_;}

  // FFT normalization factor
  const double & normFFT() const
    {return normFFT_;}

  // Dummy spectral FunctionSpace
  const atlas::FunctionSpace & spFspace() const
    {return *spFspace_;}

  // Zonal truncation
  const size_t & M() const
    {return M_;}

  // Meridional truncation
  const size_t & N() const
    {return N_;}

  // Return jk for this wavenumber
  const size_t & jk(const size_t & js) const
    {return jkVec_[js];}

  // Return jl for this wavenumber
  const size_t & jl(const size_t & js) const
    {return jlVec_[js];}

  // Return jq for this wavenumber
  const size_t & jq(const size_t & js) const
    {return jqVec_[js];}

  // Return kstar for this wavenumber
  double kstar(const size_t & js) const
    {return kstarVec_[js];}

  // Return jw for this wavenumber
  size_t jw(const size_t & js) const
    {return jwGlbVec_[js] - nwStartPerTask_[myrank_];}

  // Return spNorm for this wavenumber
  const double & spNorm(const size_t & js) const
    {return spNormVec_[js];}

  // Return spNormSumInv for this total wavenumber
  const double & spNormSumInv(const size_t & jw) const
    {return spNormSumInv_[jw];}

  // Truncation ellips
  const std::vector<size_t> & ellips() const
    {return ellips_;}

  // Global number of total wavenumber
  const size_t & nwGlb() const
    {return nwGlb_;}

  // Local number of total wavenumber
  const size_t & nw() const
    {return nw_;}

  // Local number of root total wavenumber
  const size_t & nwRoot() const
    {return nwRoot_;}

  // Vector of minimum global nw
  const std::vector<size_t> & nwStartPerTask() const
    {return nwStartPerTask_;}

  // Vector of maximum global nw
  const std::vector<size_t> & nwEndPerTask() const
    {return nwEndPerTask_;}

  // Vector of nw
  const std::vector<size_t> & nwPerTask() const
    {return nwPerTask_;}

  // Vector for root nw
  const std::vector<size_t> & nwRootPerTask() const
    {return nwRootPerTask_;}

  // Starting global nw
  const size_t & nwStart() const
    {return nwStartPerTask_[myrank_];}

  // Ending global nw
  const size_t & nwEnd() const
    {return nwEndPerTask_[myrank_];}

  // Total number of levels (sum of all levels of all active variables)
  const size_t & nvz() const
    {return nvz_;}

  // Public methods

  // Run tests
  void test(const oops::Variables &) const;

  // Spectral operations

  // Forward FFT
  void gp2sp(const atlas::FieldSet &,
             atlas::FieldSet &,
             const oops::Variables &) const;

  // Inverse FFT
  void sp2gp(const atlas::FieldSet &,
             atlas::FieldSet &,
             const oops::Variables &) const;

  // Forward FFT adjoint
  void gp2spAdj(const atlas::FieldSet &,
                    atlas::FieldSet &,
                    const oops::Variables &) const;

  // Inverse FFT adjoint
  void sp2gpAdj(const atlas::FieldSet &,
                    atlas::FieldSet &,
                    const oops::Variables &) const;

  // Create random spectral FieldSet
  void createRandomFieldSet(atlas::FieldSet &,
                            const oops::Variables &) const;

  // Apply derivative
  void derivative(const atlas::Field &,
                  atlas::Field &,
                  const std::string &,
                  const bool & adjoint = false) const;

  // Apply direct Laplacian
  void directLaplacian(atlas::Field &) const;

  // Apply inverse Laplacian
  void inverseLaplacian(atlas::Field &) const;

  // Control vector size
  size_t ctlVecSize() const
    {return ctlVecSize_;}

  // Convert control vector to spectral FieldSet
  void cv2fset(const atlas::Field &,
               atlas::FieldSet &,
               const oops::Variables &,
               const size_t &) const;

  // Convert spectral FieldSet to control vector
  void fset2cv(const atlas::FieldSet &,
               atlas::Field &,
               const oops::Variables &,
               const size_t &) const;

  // Copy spectral FieldSet
  void copyFieldSet(const atlas::FieldSet &,
                    atlas::FieldSet &,
                    const oops::Variables &) const;

  // Real kstar value
  double rkstar(const size_t &,
                const size_t &,
                const size_t &,
                const size_t &,
                const size_t &) const;


  // Integer kstar value
  double ikstar(const size_t &,
                const size_t &,
                const size_t &,
                const size_t &,
                const size_t &) const;

  // Include wavenumber in the calibration process
  bool includeWavenumber(const size_t &,
                         const size_t &) const;

  // Reduce covariances
  void reduceCov(atlas::Field &) const;

  // Scatter covariances
  void scatterCov(const std::vector<double> &,
                  atlas::Field &,
                  const bool & adjointInput = false) const;

  // Gather covariances
  void gatherCov(const atlas::Field &,
                 std::vector<double> &,
                 const bool & adjoint = false) const;

  // Filter covariances
  void filterCov(const double &,
                 atlas::Field &) const;

  // Compute covariances norm
  double normCov(const atlas::Field &) const;

  // Reduce and normalize covariance
  void reduceNormalizeCov(const size_t &,
                          atlas::Field &);

 protected:
  // Model grid geometry data
  const oops::GeometryData & gdata_;

  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Parameters
  const BifourierTransformParameters & params_;

  // FFT backend
  const std::string fftBackend_;

  // UIDs
  std::string gridUid_;
  std::string specUid_;

  // Model grid
  size_t nx_;
  size_t ny_;
  size_t nodes_;
  double dx_;
  double dy_;
  double meanLat_;

  // Spectral space
  size_t nk_;
  size_t nl_;

  // Total number of levels (sum of all levels of all active variables)
  size_t nvz_;

  // Differential operators factors
  double exwn_;
  double eywn_;

  // Normalization factor
  double normFFT_;

  // Truncation (number of positive wavenumber in y dimension for each wavenumber in x dimension)
  size_t M_;
  size_t N_;
  std::vector<size_t> ellips_;

  // Mapping and normalization
  enum Quad {
    ReRe = 0,
    ReIm = 1,
    ImRe = 2,
    ImIm = 3
  };
  struct spElem {
    size_t jk;
    size_t jl;
    Quad jq;
    double kstar;
    size_t jwGlb;
    size_t jsXDerivativeOffset;
    size_t jsYDerivativeOffset;
    size_t jt;
  };
  double jwGlbTol_;
  std::vector<spElem> spVec_;
  std::vector<size_t> spNormKL_;
  std::vector<size_t> jkVec_;
  std::vector<size_t> jlVec_;
  std::vector<size_t> jqVec_;
  std::vector<double> kstarVec_;
  std::vector<size_t> jwGlbVec_;
  std::vector<double> spNormVec_;
  std::vector<double> lapDirVec_;
  std::vector<double> lapInvVec_;
  atlas::FieldSet derivatives_;

  // Parallelization

  // Global mapping
  std::vector<size_t> sToSGlb_;

  // Grid
  size_t gridRecvSize_;
  std::vector<int> gridRecvIndex_;
  std::vector<int> gridRecvCounts_;
  std::vector<int> gridRecvDispls_;

  // Equal chunks
  size_t eqchSendSize_;
  std::vector<int> eqchSendIndex_;
  std::vector<int> eqchSendCounts_;
  std::vector<int> eqchSendDispls_;

  // Local spectral space
  std::vector<bool> truncMask_;
  std::vector<size_t> nsPerTask_;
  size_t ns_;
  size_t nsGlb_;
  std::vector<int> sCounts_;
  std::vector<int> sDispls_;
  std::vector<size_t> sMapping_;

  // Dummy spectral function space
  std::unique_ptr<atlas::FunctionSpace> spFspace_;

  // Total wavenumber
  size_t nwGlb_;
  double dwGlb_;
  size_t nw_;
  size_t nwRoot_;
  std::vector<size_t> nwStartPerTask_;
  std::vector<size_t> nwEndPerTask_;
  std::vector<size_t> nwPerTask_;
  std::vector<size_t> nwRootPerTask_;
  std::vector<int> wCounts_;
  std::vector<int> wDispls_;

  // Covariance communicators
  std::vector<bool> myJwGlb_;
  std::vector<eckit::mpi::Comm*> covRedComm_;
  std::vector<eckit::mpi::Comm*> covBcastComm_;

  // Control vector size
  size_t ctlVecSize_;

  // Spectral norm sum inverse
  std::vector<double> spNormSumInv_;

  // Private methods

  // Print
  void print(std::ostream &) const;

  // Setup global spectral space
  void setupGlobalSpectralSpace();

  // Setup parallelization
  void setupParallelizationInit();
  void setupParallelizationFinal();

  // Setup local spectral space
  void setupLocalSpectralSpace();

  // Add spectral coefficient
  void addSpectralCoefficient(const size_t &,
                              const size_t &,
                              const Quad &,
                              const size_t &,
                              const size_t &);
};

// -----------------------------------------------------------------------------

class BifourierTransformFactory;

// -----------------------------------------------------------------------------

class BifourierTransformFactory {
 public:
  static std::shared_ptr<BifourierTransformBase> create(const oops::GeometryData &,
                                                        const std::string &,
                                                        const oops::Variables &,
                                                        const BifourierTransformParameters &);

  virtual ~BifourierTransformFactory() = default;

 protected:
  explicit BifourierTransformFactory(const std::string &name);

 private:
  virtual std::shared_ptr<BifourierTransformBase> make(const oops::GeometryData &,
                                                       const std::string &,
                                                       const oops::Variables &,
                                                       const BifourierTransformParameters &) = 0;

  static std::map < std::string, BifourierTransformFactory * > & getMakers() {
    static std::map < std::string, BifourierTransformFactory * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class T>
class BifourierTransformMaker : public BifourierTransformFactory {
  std::shared_ptr<BifourierTransformBase> make(const oops::GeometryData & gdata,
                                               const std::string & gridUid,
                                               const oops::Variables & activeVars,
                                               const BifourierTransformParameters & params)
    override {
      return std::make_shared<T>(gdata, gridUid, activeVars, params);
  }

 public:
  explicit BifourierTransformMaker(const std::string & name) : BifourierTransformFactory(name) {}
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
