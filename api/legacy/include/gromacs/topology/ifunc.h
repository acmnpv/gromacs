/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#ifndef GMX_TOPOLOGY_IFUNC_H
#define GMX_TOPOLOGY_IFUNC_H

#include "gromacs/libgromacs_export.h"
#include "gromacs/math/vectypes.h"

struct t_fcdata;
struct t_graph;
union t_iparams;
struct t_mdatoms;
struct t_pbc;

/* TODO: Remove this typedef when t_ilist is removed */
typedef int t_iatom;

/* Real vector type with an additional, unused 4th element.
 * This type is used to allow aligned 4-wide SIMD loads and stores.
 */
typedef real rvec4[4];

/*
 * The function type t_ifunc() calculates one interaction, using iatoms[]
 * and iparams. Within the function the number of atoms to be used is
 * known. Within the function only the atomid part of the iatoms[] array
 * is supplied, not the type field (see also t_ilist). The function
 * returns the potential energy. If pbc==NULL the coordinates in x are
 * assumed to be such that no calculation of PBC is necessary,
 * If pbc!=NULL a full PBC calculation is performed.
 * If g!=NULL it is used for determining the shift forces.
 * With domain decomposition ddgatindex can be used for getting global
 * atom numbers for warnings and error messages.
 * ddgatindex is NULL when domain decomposition is not used.
 */

constexpr unsigned int IF_NULL       = 0;
constexpr unsigned int IF_BOND       = 1 << 0;
constexpr unsigned int IF_VSITE      = 1 << 1;
constexpr unsigned int IF_CONSTRAINT = 1 << 2;
constexpr unsigned int IF_CHEMBOND   = 1 << 3;
constexpr unsigned int IF_BTYPE      = 1 << 4;
constexpr unsigned int IF_ATYPE      = 1 << 5;
constexpr unsigned int IF_DIHEDRAL   = 1 << 6;
constexpr unsigned int IF_PAIR       = 1 << 7;
constexpr unsigned int IF_TABULATED  = 1 << 8;
constexpr unsigned int IF_LIMZERO    = 1 << 9;
/* These flags tell to some of the routines what can be done with this
 * item in the list.
 * With IF_BOND a bonded interaction will be calculated.
 * With IF_BTYPE grompp can convert the bond to a Morse potential.
 * With IF_BTYPE or IF_ATYPE the bond/angle can be converted to
 * a constraint or used for vsite parameter determination by grompp.
 * IF_LIMZERO indicates that for a bonded interaction the potential
 * does goes to zero for large distances, thus if such an interaction
 * it not assigned to any node by the domain decompostion, the simulation
 * still continue, if mdrun has been told so.
 */

struct t_interaction_function // NOLINT (clang-analyzer-optin.performance.Padding)
{
    const char* name;         /* the name of this function			*/
    const char* longname;     /* The name for printing etc.                   */
    int         nratoms;      /* nr of atoms needed for this function		*/
    int         nrfpA, nrfpB; /* number of parameters for this function.      */
                              /* this corresponds to the number of params in  */
                              /* iparams struct! (see idef.h)                 */
    /* A and B are for normal and free energy components respectively.    */
    unsigned int flags; /* Flags (see above)                            */
};

#define NRFPA(ftype) (interaction_function[(ftype)].nrfpA)
#define NRFPB(ftype) (interaction_function[(ftype)].nrfpB)
#define NRFP(ftype) (NRFPA(ftype) + NRFPB(ftype))
#define NRAL(ftype) (interaction_function[(ftype)].nratoms)

#define IS_CHEMBOND(ftype) \
    (interaction_function[(ftype)].nratoms == 2 && (interaction_function[(ftype)].flags & IF_CHEMBOND))
/* IS_CHEMBOND tells if function type ftype represents a chemical bond */

/* IS_ANGLE tells if a function type ftype represents an angle
 * Per Larsson, 2007-11-06
 */
#define IS_ANGLE(ftype) \
    (interaction_function[(ftype)].nratoms == 3 && (interaction_function[(ftype)].flags & IF_ATYPE))
#define IS_VSITE(ftype) (interaction_function[(ftype)].flags & IF_VSITE)

#define IS_TABULATED(ftype) (interaction_function[(ftype)].flags & IF_TABULATED)

/* this MUST correspond to the
   t_interaction_function[F_NRE] in src/gromacs/topology/ifunc.cpp */
enum
{
    F_BONDS,
    F_G96BONDS,
    F_MORSE,
    F_CUBICBONDS,
    F_CONNBONDS,
    F_HARMONIC,
    F_FENEBONDS,
    F_TABBONDS,
    F_TABBONDSNC,
    F_RESTRBONDS,
    F_ANGLES,
    F_G96ANGLES,
    F_RESTRANGLES,
    F_LINEAR_ANGLES,
    F_CROSS_BOND_BONDS,
    F_CROSS_BOND_ANGLES,
    F_UREY_BRADLEY,
    F_QUARTIC_ANGLES,
    F_TABANGLES,
    F_PDIHS,
    F_RBDIHS,
    F_RESTRDIHS,
    F_CBTDIHS,
    F_FOURDIHS,
    F_IDIHS,
    F_PIDIHS,
    F_TABDIHS,
    F_CMAP,
    F_GB12_NOLONGERUSED,
    F_GB13_NOLONGERUSED,
    F_GB14_NOLONGERUSED,
    F_GBPOL_NOLONGERUSED,
    F_NPSOLVATION_NOLONGERUSED,
    F_LJ14,
    F_COUL14,
    F_LJC14_Q,
    F_LJC_PAIRS_NB,
    F_LJ,
    F_BHAM,
    F_LJ_LR_NOLONGERUSED,
    F_BHAM_LR_NOLONGERUSED,
    F_DISPCORR,
    F_COUL_SR,
    F_COUL_LR_NOLONGERUSED,
    F_RF_EXCL,
    F_COUL_RECIP,
    F_LJ_RECIP,
    F_DPD,
    F_POLARIZATION,
    F_WATER_POL,
    F_THOLE_POL,
    F_ANHARM_POL,
    F_POSRES,
    F_FBPOSRES,
    F_DISRES,
    F_DISRESVIOL,
    F_ORIRES,
    F_ORIRESDEV,
    F_ANGRES,
    F_ANGRESZ,
    F_DIHRES,
    F_DIHRESVIOL,
    F_CONSTR,
    F_CONSTRNC,
    F_SETTLE,
    F_VSITE1,
    F_VSITE2,
    F_VSITE2FD,
    F_VSITE3,
    F_VSITE3FD,
    F_VSITE3FAD,
    F_VSITE3OUT,
    F_VSITE4FD,
    F_VSITE4FDN,
    F_VSITEN,
    F_COM_PULL,
    F_DENSITYFITTING,
    F_EQM,
    F_EPOT,
    F_EKIN,
    F_ETOT,
    F_ECONSERVED,
    F_TEMP,
    F_VTEMP_NOLONGERUSED,
    F_PDISPCORR,
    F_PRES,
    F_DVDL_CONSTR,
    F_DVDL,
    F_DKDL,
    F_DVDL_COUL,
    F_DVDL_VDW,
    F_DVDL_BONDED,
    F_DVDL_RESTRAINT,
    F_DVDL_TEMPERATURE, /* not calculated for now, but should just be the energy (NVT) or enthalpy (NPT), or 0 (NVE) */
    F_NRE /* This number is for the total number of energies      */
};

static inline bool IS_RESTRAINT_TYPE(int ifunc)
{
    return ifunc == F_POSRES || ifunc == F_FBPOSRES || ifunc == F_DISRES || ifunc == F_RESTRBONDS
           || ifunc == F_DISRESVIOL || ifunc == F_ORIRES || ifunc == F_ORIRESDEV
           || ifunc == F_ANGRES || ifunc == F_ANGRESZ || ifunc == F_DIHRES;
}

/* Maximum allowed number of atoms, parameters and terms in interaction_function.
 * Check kernel/toppush.c when you change these numbers.
 */
constexpr int MAXATOMLIST   = 6;
constexpr int MAXFORCEPARAM = 12;
constexpr int NR_RBDIHS     = 6;
constexpr int NR_CBTDIHS    = 6;
constexpr int NR_FOURDIHS   = 4;


constexpr t_interaction_function def_bonded(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND };
}

constexpr t_interaction_function def_pair(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_PAIR | IF_LIMZERO };
}

constexpr t_interaction_function def_bondedt(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_TABULATED };
}

constexpr t_interaction_function def_bondedtz(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_TABULATED | IF_LIMZERO };
}

constexpr t_interaction_function def_angle(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_ATYPE };
}

constexpr t_interaction_function def_dihedral(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_DIHEDRAL };
}

constexpr t_interaction_function
def_dihedral_tabulated(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_DIHEDRAL | IF_TABULATED };
}

constexpr t_interaction_function def_bond(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_CHEMBOND | IF_BTYPE };
}

constexpr t_interaction_function def_bondt(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_CHEMBOND | IF_TABULATED };
}

constexpr t_interaction_function def_bondnb(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_BOND | IF_CHEMBOND };
}

constexpr t_interaction_function def_vsite(const char* str, const char* lstr, int nra, int nrpa)
{
    return t_interaction_function{ str, lstr, nra, nrpa, 0, IF_VSITE };
}

constexpr t_interaction_function def_shk(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_CONSTRAINT };
}

constexpr t_interaction_function def_shkcb(const char* str, const char* lstr, int nra, int nrpa, int nrpb)
{
    return t_interaction_function{ str, lstr, nra, nrpa, nrpb, IF_CONSTRAINT | IF_CHEMBOND };
}

constexpr t_interaction_function def_nb(const char* str, const char* lstr, int nra, int nrpa)
{
    return t_interaction_function{ str, lstr, nra, nrpa, 0, IF_NULL };
}

constexpr t_interaction_function def_nofc(const char* str, const char* lstr)
{
    return t_interaction_function{ str, lstr, 0, 0, 0, IF_NULL };
}

const t_interaction_function interaction_function[F_NRE] = {
    def_bond("BONDS", "Bond", 2, 2, 2),
    def_bond("G96BONDS", "G96Bond", 2, 2, 2),
    def_bond("MORSE", "Morse", 2, 3, 3),
    def_bond("CUBICBONDS", "Cubic Bonds", 2, 3, 0),
    def_bondnb("CONNBONDS", "Connect Bonds", 2, 0, 0),
    def_bonded("HARMONIC", "Harmonic Pot.", 2, 2, 2),
    def_bondnb("FENEBONDS", "FENE Bonds", 2, 2, 0),
    def_bondt("TABBONDS", "Tab. Bonds", 2, 2, 2),
    def_bondedtz("TABBONDSNC", "Tab. Bonds NC", 2, 2, 2),
    def_bonded("RESTRAINTPOT", "Restraint Pot.", 2, 4, 4),
    def_angle("ANGLES", "Angle", 3, 2, 2),
    def_angle("G96ANGLES", "G96Angle", 3, 2, 2),
    def_angle("RESTRANGLES", "Restricted Angles", 3, 2, 2),
    def_angle("LINEAR_ANGLES", "Lin. Angle", 3, 2, 2),
    def_bonded("CROSS_BOND_BOND", "Bond-Cross", 3, 3, 0),
    def_bonded("CROSS_BOND_ANGLE", "BA-Cross", 3, 4, 0),
    def_angle("UREY_BRADLEY", "U-B", 3, 4, 4),
    def_angle("QANGLES", "Quartic Angles", 3, 6, 0),
    def_bondedt("TABANGLES", "Tab. Angles", 3, 2, 2),
    def_dihedral("PDIHS", "Proper Dih.", 4, 3, 3),
    def_dihedral("RBDIHS", "Ryckaert-Bell.", 4, 6, 6),
    def_dihedral("RESTRDIHS", "Restricted Dih.", 4, 2, 2),
    def_dihedral("CBTDIHS", "CBT Dih.", 4, 6, 6),
    def_dihedral("FOURDIHS", "Fourier Dih.", 4, 4, 4),
    def_dihedral("IDIHS", "Improper Dih.", 4, 2, 2),
    def_dihedral("PIDIHS", "Periodic Improper Dih.", 4, 3, 3),
    def_dihedral_tabulated("TABDIHS", "Tab. Dih.", 4, 2, 2),
    def_dihedral("CMAP", "CMAP Dih.", 5, -1, -1),
    def_nofc("GB12", "GB 1-2 Pol. (unused)"),
    def_nofc("GB13", "GB 1-3 Pol. (unused)"),
    def_nofc("GB14", "GB 1-4 Pol. (unused)"),
    def_nofc("GBPOL", "GB Polarization (unused)"),
    def_nofc("NPSOLVATION", "Nonpolar Sol. (unused)"),
    def_pair("LJ14", "LJ-14", 2, 2, 2),
    def_nofc("COUL14", "Coulomb-14"),
    def_pair("LJC14_Q", "LJC-14 q", 2, 5, 0),
    def_pair("LJC_NB", "LJC Pairs NB", 2, 4, 0),
    def_nb("LJ_SR", "LJ (SR)", 2, 2),
    def_nb("BHAM", "Buck.ham (SR)", 2, 3),
    def_nofc("LJ_LR", "LJ (unused)"),
    def_nofc("BHAM_LR", "B.ham (unused)"),
    def_nofc("DISPCORR", "Disper. corr."),
    def_nofc("COUL_SR", "Coulomb (SR)"),
    def_nofc("COUL_LR", "Coul (unused)"),
    def_nofc("RF_EXCL", "RF excl."),
    def_nofc("COUL_RECIP", "Coul. recip."),
    def_nofc("LJ_RECIP", "LJ recip."),
    def_nofc("DPD", "DPD"),
    def_bondnb("POLARIZATION", "Polarization", 2, 1, 0),
    def_bonded("WATERPOL", "Water Pol.", 5, 6, 0),
    def_bonded("THOLE", "Thole Pol.", 4, 3, 0),
    def_bondnb("ANHARM_POL", "Anharm. Pol.", 2, 3, 0),
    def_bonded("POSRES", "Position Rest.", 1, 3, 3),
    def_bonded("FBPOSRES", "Flat-bottom posres", 1, 3, 0),
    def_bonded("DISRES", "Dis. Rest.", 2, 6, 0),
    def_nofc("DISRESVIOL", "D.R.Viol. (nm)"),
    def_bonded("ORIRES", "Orient. Rest.", 2, 6, 0),
    def_nofc("ORDEV", "Ori. R. RMSD"),
    def_bonded("ANGRES", "Angle Rest.", 4, 3, 3),
    def_bonded("ANGRESZ", "Angle Rest. Z", 2, 3, 3),
    def_bonded("DIHRES", "Dih. Rest.", 4, 3, 3),
    def_nofc("DIHRESVIOL", "Dih. Rest. Viol."), /* obsolete */
    def_shkcb("CONSTR", "Constraint", 2, 1, 1),
    def_shk("CONSTRNC", "Constr. No Conn.", 2, 1, 1),
    def_shkcb("SETTLE", "Settle", 3, 2, 0),
    def_vsite("VSITE1", "Virtual site 1", 2, 0),
    def_vsite("VSITE2", "Virtual site 2", 3, 1),
    def_vsite("VSITE2FD", "Virtual site 2fd", 3, 1),
    def_vsite("VSITE3", "Virtual site 3", 4, 2),
    def_vsite("VSITE3FD", "Virtual site 3fd", 4, 2),
    def_vsite("VSITE3FAD", "Virtual site 3fad", 4, 2),
    def_vsite("VSITE3OUT", "Virtual site 3out", 4, 3),
    def_vsite("VSITE4FD", "Virtual site 4fd", 5, 3),
    def_vsite("VSITE4FDN", "Virtual site 4fdn", 5, 3),
    def_vsite("VSITEN", "Virtual site N", 2, 2),
    def_nofc("COM_PULL", "COM Pull En."),
    def_nofc("DENSITYFIT", "Density fitting"),
    def_nofc("EQM", "Quantum En."),
    def_nofc("EPOT", "Potential"),
    def_nofc("EKIN", "Kinetic En."),
    def_nofc("ETOT", "Total Energy"),
    def_nofc("ECONS", "Conserved En."),
    def_nofc("TEMP", "Temperature"),
    def_nofc("VTEMP", "Vir. Temp. (not used)"),
    /* Note that pressure names can not be more than 8 char's,
     * because " (bar)" is appended to them.
     */
    def_nofc("PDISPCORR", "Pres. DC"),
    def_nofc("PRES", "Pressure"),
    def_nofc("DH/DL_CON", "dH/dl constr."), /* obsolete */
    def_nofc("DV/DL", "dVremain/dl"),
    def_nofc("DK/DL", "dEkin/dl"),
    def_nofc("DVC/DL", "dVcoul/dl"),
    def_nofc("DVV/DL", "dVvdw/dl"),
    def_nofc("DVB/DL", "dVbonded/dl"),
    def_nofc("DVR/DL", "dVrestraint/dl"),
    def_nofc("DVT/DL", "dVtemperature/dl")
};

#endif
