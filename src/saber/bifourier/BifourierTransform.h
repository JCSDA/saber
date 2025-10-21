/*
 * (C) Copyright 2025 Meteorologisk Institutt
 *
 */

#pragma once

#include <fftw3.h>

#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "eckit/mpi/Comm.h"
#include "oops/util/Printable.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"

namespace saber {
namespace bifourier {

// -----------------------------------------------------------------------------

class BifourierTransform : public util::Printable {
 public:
  static const std::string classname()
    {return "saber::bifourier::BifourierTransform";}

  // Constructors

  // Create spectral transform
  BifourierTransform(const oops::GeometryData &,
                     const std::string &,
                     const oops::Variables &,
                     const eckit::Configuration &);

  // Destructor
  ~BifourierTransform();

  // Accessors

  // Geometry data
  const oops::GeometryData & geometryData() const
    {return gdata_;}

  // UIDs
  const std::string & gridUid() const
    {return gridUid_;}
  const std::string & specUid() const
    {return specUid_;}

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

  // Zonal wavenumbers size
  const size_t & nk() const
    {return nk_;}

  // Meridional wavenumbers size
  const size_t & nl() const
    {return nl_;}

  // Return jk for this wavenumber
  const size_t & jk(const size_t & js) const
    {return jkVec_[js];}

  // Return jl for this wavenumber
  const size_t & jl(const size_t & js) const
    {return jlVec_[js];}

  // Return jq for this wavenumber
  const size_t & jq(const size_t & js) const
    {return jqVec_[js];}

  // Return jw for this wavenumber
  const size_t jw(const size_t & js) const
    {return jwGlbVec_[js] - nwStartPerTask_[myrank_];}

  // Return spNorm for this wavenumber
  const double & spNorm(const size_t & js) const
    {return spNormVec_[js];}

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

  // Communication vectors
  const std::vector<int> & wCounts() const
    {return wCounts_;}
  const std::vector<int> & wDispls() const
    {return wDispls_;}

  // Total number of levels (sum of all levels of all active variables)
  const size_t & nvz() const
    {return nvz_;}

  // Spectral operations

  // Forward FFT
  void gp2sp(const atlas::FieldSet &,
             atlas::FieldSet &,
             const oops::Variables &) const;

  // Inverse FFT
  void sp2gp(const atlas::FieldSet &,
             atlas::FieldSet &,
             const oops::Variables &) const;

  // Create random spectral FieldSet
  void createRandomSpectralFieldSet(atlas::FieldSet &,
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

  // kstar value
  double kstar(const size_t &,
               const size_t &,
               const size_t &,
               const size_t &,
               const size_t &) const;

 private:
  // Model grid geometry data
  const oops::GeometryData & gdata_;

  // Communicator
  const eckit::mpi::Comm & comm_;
  size_t myrank_;

  // Parameters
  const eckit::Configuration & params_;

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
  std::vector<size_t> jwGlbVec_;
  std::vector<double> spNormVec_;
  std::vector<double> lapDirVec_;
  std::vector<double> lapInvVec_;
  atlas::FieldSet derivatives_;

  // Parallelization
  std::vector<size_t> nyPerTask_;
  std::vector<size_t> nkPerTask_;
  std::vector<size_t> nsPerTask_;
  std::vector<size_t> sToSGlb_;

  // Rows <=> grid
  size_t rowsSendSize_;
  size_t gridRecvSize_;
  std::vector<int> gridRecvIndex_;
  std::vector<int> rowsSendIndex_;
  std::vector<int> rowsSendCounts_;
  std::vector<int> rowsSendDispls_;
  std::vector<int> gridRecvCounts_;
  std::vector<int> gridRecvDispls_;

  // Columns <=> rows
  size_t colsSendSize_;
  size_t rowsRecvSize_;
  std::vector<int> rowsRecvIndex_;
  std::vector<int> colsSendIndex_;
  std::vector<int> colsSendCounts_;
  std::vector<int> colsSendDispls_;
  std::vector<int> rowsRecvCounts_;
  std::vector<int> rowsRecvDispls_;

  // Equal chunks <=> columns
  size_t eqchSendSize_;
  size_t colsRecvSize_;
  std::vector<int> colsRecvIndex_;
  std::vector<std::vector<int>> colsJq_;
  std::vector<int> eqchSendIndex_;
  std::vector<int> eqchSendCounts_;
  std::vector<int> eqchSendDispls_;
  std::vector<int> colsRecvCounts_;
  std::vector<int> colsRecvDispls_;

  // Local spectral space
  std::vector<bool> truncMask_;
  size_t ns_;
  size_t nsGlb_;
  std::vector<int> sCounts_;
  std::vector<int> sDispls_;
  std::vector<size_t> sMapping_;

  // Dummy spectral function space
  std::unique_ptr<atlas::FunctionSpace> spFspace_;

  // Control vector size
  size_t ctlVecSize_;

  // Total wavenumber
  size_t nwGlb_;
  size_t nw_;
  size_t nwRoot_;
  std::vector<size_t> nwStartPerTask_;
  std::vector<size_t> nwEndPerTask_;
  std::vector<size_t> nwPerTask_;
  std::vector<size_t> nwRootPerTask_;
  std::vector<int> wCounts_;
  std::vector<int> wDispls_;

  // Rows FFT
  fftw_plan rowsPlan_r2c_;
  fftw_plan rowsPlan_c2r_;
  double *rowsBufR_ = nullptr;
  fftw_complex *rowsBufC_ = nullptr;

  // Columns FFT
  fftw_plan colsPlan_r2c_;
  fftw_plan colsPlan_c2r_;
  double *colsBufR_ = nullptr;
  fftw_complex *colsBufC_ = nullptr;

  // Private methods

  // Print
  void print(std::ostream &) const;

  // Setup global spectral space
  void setupGlobalSpectralSpace();

  // Setup parallelization
  void setupParallelization();

  // Setup local spectral space
  void setupLocalSpectralSpace();

  // Setup FFT
  void setupFFT();
  void cleanupFFT();

  // Add spectral coefficient
  void addSpectralCoefficient(const size_t &,
                              const size_t &,
                              const Quad &,
                              const size_t &,
                              const size_t &);
};

// -----------------------------------------------------------------------------

}  // namespace bifourier
}  // namespace saber
