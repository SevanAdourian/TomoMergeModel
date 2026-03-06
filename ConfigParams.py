# imports here


class ConfigParams:
    def __init__(
        self,
        reg_lmax: int,
        win_lmax: int,
        win_eff_lmax: int,
        lon_bounds: tuple[float],
        lat_bounds: tuple[float],
        win_type: str,
    ):
        """Initialize configuration parameters for merging regional and global models."""

        self.reg_lmax = reg_lmax
        self.win_lmax = win_lmax
        self.win_eff_lmax = win_eff_lmax
        self.lon_min_mask, self.lon_max_mask = lon_bounds
        self.lat_min_mask, self.lat_max_mask = lat_bounds
        self.win_type = win_type
