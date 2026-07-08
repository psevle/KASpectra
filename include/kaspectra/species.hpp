#pragma once

// Shared particle species enum for both pp and pgamma
// production channels
// e+/e- and nu/nubar kept as 7 distinct entries because pgamma needs
// them distinct
// For pgamma, there are distinct fitted tables for particle/antiparticle
// whereas pp considers them the same

namespace kaspectra {

    enum class Species {
        Gamma,
        Positron,
        Electron,
        NuMu,
        NuMuBar,
        NuE,
        NuEBar
    };

} // namespace kaspectra