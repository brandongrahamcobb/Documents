---

# Discord Bot Project

This repository contains a highly customizable Discord bot that offers a wide range of functionalities, including molecule analysis, information retrieval, and user interaction enhancements. Below is an overview of the bot's main features, commands, and structure.

## Features

### Analog Command
The `analog` command allows users to upload a list of molecular structures, separated by spaces. The bot initiates an active dialog between Lucy (the bot's character) and the user to calculate the Tanimoto Similarity of the molecules in the list. This feature helps identify which molecule is most structurally similar within the list, and it can deduce potential interactions between molecules based on their structural similarity. This command is currently a regular command but will soon be moved to a slash command.

### Draw Command
The `draw` command is a hybrid command that allows users to visualize chemical structures. It integrates with RDKit and other chemistry-related tools to generate images of the molecular structures.

### Slash Commands
All user functions, except for `analog` and `draw`, have been migrated to slash commands. This makes the bot more interactive and modern, utilizing Discord's latest command features for a smoother user experience.

### Helpers File
The bot's `helpers.py` file acts as an extension to the core of the bot, providing a wide array of functionalities, including:

- **Bible Script Recall**: Retrieve specific Bible verses on command.
- **Wikipedia Searching**: Search Wikipedia for articles and retrieve relevant information, including the first image on the page.
- **Google Searches**: Perform Google searches and return the top results directly within Discord.
- **Google Translate**: Translate text between different languages using Google Translate API.
- **Molecule Analysis**: Analyze molecular structures, calculate Tanimoto Similarity, and perform various chemical operations.
- **Tag System**: A simplified tag system to store and retrieve user-generated content.

### Additional Features
- **Wikipedia Plant Category**: Retrieve random Wikipedia pages categorized under plants and display relevant information.
- **Cooldown and Heatup Management**: Dynamic command cooldowns to manage how often users can invoke certain commands.
- **Image Processing**: Manipulate images in bytes using the Python Imaging Library (PIL) within a Discord bot cog.
- **Database Operations**: Integration with `aiomysql` for efficient database operations.
- **Dynamic Bot Control**: Start, stop, and restart the bot dynamically, including server icon management.

## Project Structure

- **cogs/**: Contains all the bot's cogs, including AdminCog, GameCog, AICog, and UserCog.
- **utils/**: Contains utility scripts, including the `helpers.py` file with core extensions.
- **main.py**: The entry point of the bot, initializing all cogs and starting the bot.

## Installation

To install and run the bot, follow the steps below:

1. Download the latest wheel release from the [Documents repository](https://github.com/brandongrahamcobb/Documents/releases).

2. Install the wheel using pip:
   ```bash
   pip install path/to/your/CobbBrandonGraham_lucy_010924-1.0.0-py3-none-any.whl
   ```
3. Run the bot:
   ```bash
   python - bot.main
   ```

## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the GPL License.

---
