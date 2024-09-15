    @commands.Cog.listener()
    async def on_command_error(self, ctx: commands.Context, error):
        if isinstance(error, commands.CommandOnCooldown):
            await ctx.send(f'Command is on cooldown. Try again in {error.retry_after:.2f} seconds.')
        else:
            raise error

