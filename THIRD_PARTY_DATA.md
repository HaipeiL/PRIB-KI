# Third-party data and attribution

## TargetTrack

PRIB's optional TargetTrack workflow uses the Protein Structure Initiative –
TargetTrack 2000–2017 archive as a historical descriptive source.

- Source: https://doi.org/10.5281/zenodo.821654
- Archive: TargetTrack-1Jul2017.tar.gz
- Expected MD5: 200012a8a2a11ffd7e370ed142df36c3
- License: Creative Commons Attribution-ShareAlike 4.0 International
  (CC BY-SA 4.0), https://creativecommons.org/licenses/by-sa/4.0/
- Creators: Helen M. Berman, Margaret J. Gabanyi, Andrei Kouranov,
  David I. Micallef, John Westbrook, and the Protein Structure Initiative
  network of investigators.

Suggested attribution:

> This analysis uses the Protein Structure Initiative TargetTrack 2000–2017
> archive created by Helen M. Berman, Margaret J. Gabanyi, Andrei Kouranov,
> David I. Micallef, John Westbrook, and the Protein Structure Initiative
> network of investigators, available at https://doi.org/10.5281/zenodo.821654
> under CC BY-SA 4.0. PRIB transformed and normalized source records; the
> original creators do not endorse PRIB or this analysis.

### PRIB transformations

The repository contains code and synthetic fixtures only. The TargetTrack
workflow is designed to perform streaming XML inspection, field normalization,
canonical stage mapping, deterministic stop classification, trial and target
aggregation, and aggregate progression reporting. It excludes personal contact
fields from public outputs.

The 832.9 MB source archive and any row-level normalized TargetTrack database
are not included in this repository. A public release containing all or a
substantial part of a TargetTrack-derived database must retain attribution and
meet the source licence's ShareAlike requirements. Code remains separate from
data licensing. This notice is an engineering control, not legal advice.
