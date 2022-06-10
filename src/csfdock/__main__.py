"""Command-line interface."""
import click


@click.command()
@click.version_option()
def main() :
    """CsfDock."""


if __name__ == "__main__":
    main(prog_name="csfdock")  # pragma: no cover
