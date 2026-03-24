# imports here


class ConfigParams:
    def __init__(
        self,
        reg_lmax: int,
        win_lmax: int,
        lon_bounds: tuple[float, float],
        lat_bounds: tuple[float, float],
        win_type: str,
        mask_mode: str = "bounds",
        mask_continents: list[str] | tuple[str, ...] | None = None,
        mask_target: str = "land",
        mask_resolution: str = "110m",
    ):
        """Initialize configuration parameters for merging regional and global models."""

        if reg_lmax >= 0:
            self.reg_lmax = reg_lmax
        else:
            raise ValueError("reg_lmax must be non-negative")
        if win_lmax >= 0:
            self.win_lmax = win_lmax
        else:
            raise ValueError("win_lmax must be non-negative")
        if mask_mode not in ["bounds", "continent"]:
            raise ValueError("mask_mode must be 'bounds' or 'continent'")
        self.mask_mode = mask_mode

        if self.mask_mode == "bounds":
            if len(lon_bounds) != 2:
                raise ValueError("lon_bounds should only have 2 values")
            elif lon_bounds[0] > lon_bounds[1]:
                raise ValueError("lon_bounds[0] can't be greater than lon_bounds[1]")
            else:
                for bound in lon_bounds:
                    if not (0 <= bound <= 360):
                        raise ValueError(f"lon_bounds value {bound} is not between 0 and 360")
                self.lon_min_mask, self.lon_max_mask = lon_bounds
            if len(lat_bounds) != 2:
                raise ValueError("lat_bounds should only have 2 values")
            elif lat_bounds[0] > lat_bounds[1]:
                raise ValueError("lat_bounds[0] can't be greater than lat_bounds[1]")
            else:
                for bound in lat_bounds:
                    if not (-90 <= bound <= 90):
                        raise ValueError(f"lat_bounds value {bound} is not between -90 and 90")
                self.lat_min_mask, self.lat_max_mask = lat_bounds
        else:
            self.lon_min_mask = None
            self.lon_max_mask = None
            self.lat_min_mask = None
            self.lat_max_mask = None

        if mask_target not in ["land", "ocean"]:
            raise ValueError("mask_target must be 'land' or 'ocean'")
        self.mask_target = mask_target
        self.mask_resolution = mask_resolution
        self.mask_continents = list(mask_continents) if mask_continents is not None else []
        
        if win_type not in ['spherical', 'rectangular']:
            raise ValueError("win_type must be 'spherical' or 'rectangular'")
        else:
            self.win_type = win_type
