from discord.ext import commands

class AICog(commands.Cog):
    def __init__(self, bot):
        self.bot = bot

async def setup(bot):
    await bot.add_cog(AICog(bot))
