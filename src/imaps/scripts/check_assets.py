"""Check that iMaps has all the assets that it needs to operate."""
import resdk
from imaps.base.constants import IMAPS_URL
from imaps.base.constants.assets import (
    ANNOTATION,
    GENOME,
    MASKED_GENOME,
    SEGMENT,
    SPECIES,
    STAR_INDEX,
    TRNA_RRNA_INDEX,
    TRNA_RRNA_SEQ,
)


def check_asset(res, asset_data, errors):
    """Check single asset."""
    # 1. Check that object exists:
    if not asset_data:
        return

    try:
        obj = res.data.get(asset_data["slug"])
    except LookupError:
        errors.append("Object with slug {} does not exist.".format(asset_data["slug"]))
        return

    # 2. Check that object has correct name:
    if obj.name != asset_data["name"]:
        errors.append(
            'Object with slug {} has wrong name: "{}" != "{}".'.format(obj.slug, obj.name, asset_data["name"])
        )


def main():
    """Run."""
    res = resdk.Resolwe(url=IMAPS_URL)
    res.login()

    errors = []

    for species in SPECIES:
        check_asset(res, GENOME[species], errors)
        check_asset(res, MASKED_GENOME.get(species, ""), errors)
        check_asset(res, ANNOTATION[species], errors)
        check_asset(res, SEGMENT[species], errors)
        check_asset(res, STAR_INDEX[species], errors)
        check_asset(res, TRNA_RRNA_SEQ[species], errors)
        check_asset(res, TRNA_RRNA_INDEX[species], errors)

    if errors:
        for err in errors:
            print(err)
        raise ValueError("See errors above.")

    print("All good, assets as expected.")


if __name__ == "__main__":
    main()
